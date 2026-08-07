"""BIGSdb REST API client for scheme retrieval.

This module provides a client for accessing BIGSdb REST endpoints to retrieve
scheme metadata, loci information, allele sequences, and profile tables.
It handles automatic pagination, temporal filtering, and captures provenance
timestamps for all responses.
"""

from dataclasses import dataclass
from datetime import date, datetime, timezone
from pathlib import Path
from typing import List, Optional, Dict, Any, Tuple, Union
import csv
import io
import re
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from torchbase.conversions.log import get_logger, TRACE

_log = get_logger("bigsdb")


# Custom exceptions


class BIGSdbError(Exception):
    """Base exception for BIGSdb client errors."""

    pass


class BIGSdbNetworkError(BIGSdbError):
    """Exception for network-related errors."""

    pass


class BIGSdbValidationError(BIGSdbError):
    """Exception for data validation errors."""

    pass


# Dataclasses for structured responses


@dataclass
class DatabaseInfo:
    """A BIGSdb sequence-definition database.

    Attributes:
        name: Database identifier used in API paths (e.g., pubmlst_neisseria_seqdef)
        description: Human-readable description
    """

    name: str
    description: str


@dataclass
class SchemeInfo:
    """A scheme entry returned by the /schemes listing endpoint.

    Attributes:
        scheme_id: Numeric identifier
        description: Human-readable scheme name
    """

    scheme_id: int
    description: str


@dataclass
class SchemeMetadata:
    """Metadata about a typing scheme.

    Attributes:
        scheme_id: Numeric identifier of the scheme
        name: Name/description of the scheme
        description: Full description text
        last_updated: Timestamp when scheme was last updated
    """

    scheme_id: int
    name: str
    description: Optional[str] = None
    last_updated: Optional[datetime] = None


@dataclass
class LocusData:
    """Information about a single locus in a scheme.

    Attributes:
        locus_id: Identifier string for the locus (e.g., "adk")
        locus_name: Human-readable name of the locus
        alleles_count: Number of known alleles for this locus
        last_updated: Timestamp when locus data was last updated
    """

    locus_id: str
    locus_name: str
    alleles_count: int
    last_updated: Optional[datetime] = None


@dataclass
class ProfileTable:
    """Allelic profile table data.

    Attributes:
        profiles: List of profile records as dictionaries
        row_count: Number of profiles in the table
        last_updated: Timestamp when profiles were last updated
    """

    profiles: List[Dict[str, str]]
    row_count: int
    last_updated: Optional[datetime] = None


@dataclass
class SchemeData:
    """Complete scheme data including metadata, loci, and profiles.

    Attributes:
        metadata: SchemeMetadata with scheme information
        loci: List of LocusData for each locus in scheme
        profiles: ProfileTable with allelic profiles
        cutoff_date: Upper date bound used when fetching allele sequences
    """

    metadata: SchemeMetadata
    loci: List[LocusData]
    profiles: Optional[ProfileTable] = None
    cutoff_date: Optional[date] = None


class BIGSdbClient:
    """Client for BIGSdb REST API.

    Handles communication with BIGSdb REST endpoints to retrieve typing scheme
    information. Automatically handles pagination, temporal filtering, and
    captures last_updated timestamps for provenance tracking.

    Attributes:
        base_url: Base URL for BIGSdb REST API
        timeout: Request timeout in seconds (default: 30)
    """

    def __init__(
        self,
        base_url: str,
        timeout: int = 30,
        verify: Union[bool, str] = True,
        retries: int = 5,
    ):
        """Initialize BIGSdb client.

        Args:
            base_url: Base URL of BIGSdb REST API
            timeout: Read timeout in seconds (connect timeout is always 10s)
            verify: SSL verification — True (default), False (dangerous, for broken VPN
                environments), or a path string to a CA bundle file
            retries: Number of retry attempts for transient failures
        """
        self.base_url = base_url
        self.timeout = timeout
        self.verify = verify
        self._session = self._make_session(retries)

    def _make_session(self, retries: int) -> requests.Session:
        session = requests.Session()
        retry_strategy = Retry(
            total=retries,
            backoff_factor=1.5,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["GET"],
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("https://", adapter)
        session.mount("http://", adapter)
        return session

    def _parse_datetime(self, datetime_str: str) -> datetime:
        """Parse ISO format datetime string to datetime object.

        Args:
            datetime_str: ISO format datetime string

        Returns:
            Parsed datetime object

        Raises:
            BIGSdbValidationError: If datetime string cannot be
                parsed
        """
        try:
            iso_str = datetime_str.replace("Z", "+00:00")
            return datetime.fromisoformat(iso_str)
        except (ValueError, AttributeError) as e:
            msg = f"Cannot parse datetime: {datetime_str}"
            raise BIGSdbValidationError(msg) from e

    def _make_request(
        self,
        method: str,
        endpoint: str,
        params: Optional[Dict[str, Any]] = None,
        expect_json: bool = True,
    ) -> Any:
        """Make HTTP request to BIGSdb API.

        Args:
            method: HTTP method (GET, POST, etc.)
            endpoint: API endpoint path
            params: Query parameters
            expect_json: Whether response should be parsed as JSON

        Returns:
            Response data (dict if JSON, str if text)

        Raises:
            BIGSdbNetworkError: For network-related errors
            BIGSdbError: For HTTP errors or API errors
            BIGSdbValidationError: For response parsing errors
        """
        url = f"{self.base_url}/{endpoint}"

        _log.log(TRACE, "HTTP GET %s (params=%s)", url, params or {})
        try:
            response = self._session.get(
                url, params=params, timeout=(10, self.timeout), verify=self.verify
            )
        except (
            requests.ConnectionError,
            requests.Timeout,
            requests.RequestException,
            ConnectionError,
            TimeoutError,
        ) as e:
            raise BIGSdbNetworkError(f"Network error: {str(e)}") from e

        if _log.isEnabledFor(TRACE):
            try:
                content_len = len(response.content)
            except TypeError:
                content_len = -1
            _log.log(TRACE, "  → HTTP %s, %d bytes", response.status_code, content_len)
        if response.status_code >= 400:
            msg = f"API error {response.status_code}: {response.reason}"
            raise BIGSdbError(msg)

        if expect_json:
            try:
                return response.json()
            except ValueError as e:
                raise BIGSdbValidationError(
                    f"Invalid JSON response: {str(e)}"
                ) from e
        else:
            return response.text

    def _fetch_scheme_metadata(
        self, database: str, scheme_id: int
    ) -> SchemeMetadata:
        """Fetch metadata for a scheme.

        Args:
            database: Database identifier (e.g., "pubmlst")
            scheme_id: Scheme ID number

        Returns:
            SchemeMetadata object

        Raises:
            BIGSdbError: If scheme cannot be fetched
        """
        endpoint = f"db/{database}/schemes/{scheme_id}"
        data = self._make_request("GET", endpoint)

        # API returns a flat dict, not {"records": [...]}
        if not data.get("id") and not data.get("description"):
            msg = f"Scheme {scheme_id} not found in database {database}"
            raise BIGSdbError(msg)

        last_updated = None
        if data.get("last_updated"):
            try:
                last_updated = datetime.strptime(
                    data["last_updated"][:10], "%Y-%m-%d"
                ).replace(tzinfo=timezone.utc)
            except ValueError:
                pass

        return SchemeMetadata(
            scheme_id=int(data.get("id", scheme_id)),
            name=data.get("description", f"Scheme {scheme_id}"),
            description=data.get("description"),
            last_updated=last_updated,
        )

    def _fetch_loci(
        self,
        database: str,
        scheme_id: int,
        added_after: Optional[datetime] = None,
        updated_after: Optional[datetime] = None,
    ) -> List[LocusData]:
        """Fetch loci for a scheme with automatic pagination.

        Args:
            database: Database identifier
            scheme_id: Scheme ID number
            added_after: Filter loci added after this datetime
            updated_after: Filter loci updated after this datetime

        Returns:
            List of LocusData objects

        Raises:
            BIGSdbError: If loci cannot be fetched
        """
        loci = []
        page = 1
        params = {}

        if added_after:
            params["added_after"] = added_after.isoformat()
        if updated_after:
            params["updated_after"] = updated_after.isoformat()

        # Loci endpoint returns all loci in one response (no pagination in practice).
        # Response: {"loci": ["https://.../loci/abcZ", ...], "records": N}
        endpoint = f"db/{database}/schemes/{scheme_id}/loci"
        _log.debug("Fetching loci for scheme %s", scheme_id)
        data = self._make_request("GET", endpoint, params=params if params else None)

        for locus_url in data.get("loci", []):
            locus_id = locus_url.rstrip("/").split("/")[-1]
            loci.append(LocusData(
                locus_id=locus_id,
                locus_name=locus_id,
                alleles_count=0,
            ))

        _log.debug("  scheme %s: %d loci", scheme_id, len(loci))
        return loci

    def _fetch_profiles(
        self,
        database: str,
        scheme_id: int,
        added_after: Optional[datetime] = None,
        updated_after: Optional[datetime] = None,
    ) -> ProfileTable:
        """Fetch allelic profiles for a scheme as CSV.

        Args:
            database: Database identifier
            scheme_id: Scheme ID number
            added_after: Filter profiles added after this datetime
            updated_after: Filter profiles updated after this datetime

        Returns:
            ProfileTable object

        Raises:
            BIGSdbError: If profiles cannot be fetched
        """
        params = {}

        if added_after:
            params["added_after"] = added_after.isoformat()
        if updated_after:
            params["updated_after"] = updated_after.isoformat()

        _log.debug("Fetching profiles for scheme %s", scheme_id)
        endpoint = f"db/{database}/schemes/{scheme_id}/profiles_csv"
        csv_text = self._make_request(
            "GET", endpoint, params=params, expect_json=False
        )

        # Parse CSV content
        profiles = []
        reader = csv.DictReader(io.StringIO(csv_text), delimiter="\t")
        if reader is None or reader.fieldnames is None:
            msg = "Invalid CSV response"
            raise BIGSdbValidationError(msg)

        for row in reader:
            profiles.append(row)

        _log.debug("  profiles: %d rows", len(profiles))
        # Capture current timestamp for provenance
        last_updated = datetime.now(timezone.utc)

        return ProfileTable(
            profiles=profiles,
            row_count=len(profiles),
            last_updated=last_updated,
        )

    def _fetch_alleles_fasta(
        self,
        database: str,
        locus_id: str,
        dest_path: Path,
        cutoff_date: Optional[date] = None,
    ) -> Path:
        """Fetch allele sequences for a locus and stream them to dest_path.

        Uses the alleles_fasta endpoint, which for unauthenticated requests
        already restricts results to alleles submitted on or before 2024-12-31
        (PubMLST's licensing cutoff).  The cutoff_date parameter is accepted
        for API compatibility but the server enforces the restriction itself.

        Args:
            database: Database identifier (e.g., "pubmlst_neisseria_seqdef")
            locus_id: Locus identifier string
            dest_path: Destination file path; parent directories are created
            cutoff_date: Accepted for compatibility; server enforces this limit

        Returns:
            dest_path

        Raises:
            BIGSdbNetworkError: For network-related errors
            BIGSdbError: If alleles cannot be fetched
        """
        dest_path = Path(dest_path)
        dest_path.parent.mkdir(parents=True, exist_ok=True)

        url = f"{self.base_url}/db/{database}/loci/{locus_id}/alleles_fasta"
        _log.log(TRACE, "HTTP GET %s (stream)", url)
        try:
            resp = self._session.get(
                url, stream=True, timeout=(10, self.timeout), verify=self.verify
            )
        except (
            requests.ConnectionError,
            requests.Timeout,
            requests.RequestException,
            ConnectionError,
            TimeoutError,
        ) as e:
            raise BIGSdbNetworkError(f"Network error fetching {locus_id}: {e}") from e

        if resp.status_code >= 400:
            raise BIGSdbError(
                f"API error {resp.status_code} fetching {locus_id}: {resp.reason}"
            )

        total = int(resp.headers.get("content-length", 0))

        try:
            from tqdm import tqdm
            bar: Any = tqdm(
                total=total or None,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                desc=locus_id,
                leave=False,
            )
        except ImportError:
            bar = None

        allele_count = 0
        try:
            with open(dest_path, "wb") as f:
                for chunk in resp.iter_content(chunk_size=1 << 16):
                    allele_count += chunk.count(b">")
                    f.write(chunk)
                    if bar is not None:
                        bar.update(len(chunk))
        finally:
            if bar is not None:
                bar.close()

        _log.debug("  %s: %d alleles → %s", locus_id, allele_count, dest_path)
        _log.log(TRACE, "  dest: %s (%d bytes)", dest_path, dest_path.stat().st_size)
        return dest_path

    def list_databases(self) -> List[DatabaseInfo]:
        """List all databases exposed by this BIGSdb instance.

        Hits the API root and returns every entry that looks like a
        sequence-definition database (name ends with '_seqdef').  Isolates
        databases are skipped because they do not carry allele sequences or
        profiles.

        Returns:
            List of DatabaseInfo objects for seqdef databases

        Raises:
            BIGSdbError: If the database list cannot be fetched
        """
        _log.info("Enumerating databases from %s", self.base_url)
        # GET /db returns a list of organism groups, each with a nested
        # "databases" array: [{"name": "...", "databases": [{...}, ...]}, ...]
        data = self._make_request("GET", "db")

        databases = []
        groups = data if isinstance(data, list) else data.get("databases", [])
        for group in groups:
            for entry in group.get("databases", [group] if "name" in group else []):
                name = entry.get("name", "")
                if name.endswith("_seqdef"):
                    databases.append(
                        DatabaseInfo(name=name, description=entry.get("description", name))
                    )

        _log.debug("  found %d seqdef databases", len(databases))
        for db in databases:
            _log.log(TRACE, "    database: %s", db.name)
        return databases

    def list_schemes(self, database: str) -> List[SchemeInfo]:
        """List all schemes in a database.

        Args:
            database: Database identifier (e.g., pubmlst_neisseria_seqdef)

        Returns:
            List of SchemeInfo objects

        Raises:
            BIGSdbError: If the scheme list cannot be fetched
        """
        _log.debug("Listing schemes in %s", database)
        endpoint = f"db/{database}/schemes"
        data = self._make_request("GET", endpoint)

        raw = data.get("schemes", data.get("records", []))
        schemes = []
        for entry in raw:
            if isinstance(entry, str):
                # URL like ".../schemes/1"
                m = re.search(r"/schemes/(\d+)", entry)
                if m:
                    scheme_id = int(m.group(1))
                    schemes.append(SchemeInfo(scheme_id=scheme_id, description=f"Scheme {scheme_id}"))
            else:
                href = entry.get("scheme", entry.get("href", ""))
                m = re.search(r"/schemes/(\d+)", href)
                if not m:
                    continue
                scheme_id = int(m.group(1))
                description = entry.get("description", entry.get("name", f"Scheme {scheme_id}"))
                schemes.append(SchemeInfo(scheme_id=scheme_id, description=description))

        for s in schemes:
            _log.log(TRACE, "    scheme %d: %s", s.scheme_id, s.description)
        return schemes

    def fetch_scheme(
        self,
        database: str,
        scheme_id: int,
        added_after: Optional[datetime] = None,
        updated_after: Optional[datetime] = None,
        cutoff_date: Optional[date] = None,
    ) -> SchemeData:
        """Fetch complete scheme data.

        Retrieves scheme metadata, loci, and profiles in one call.
        Automatically handles pagination for loci and applies
        temporal filters.

        Args:
            database: Database identifier (e.g., "pubmlst")
            scheme_id: Numeric scheme identifier
            added_after: Filter for records added after this datetime
            updated_after: Filter for records updated after this datetime
            cutoff_date: Exclude allele sequences entered after this date;
                use date(2024, 12, 31) to stay within freely-licensed data

        Returns:
            SchemeData object containing all scheme information

        Raises:
            BIGSdbError: If scheme cannot be fetched
            BIGSdbNetworkError: For network errors
            BIGSdbValidationError: For data validation errors
        """
        metadata = self._fetch_scheme_metadata(database, scheme_id)
        loci = self._fetch_loci(
            database,
            scheme_id,
            added_after=added_after,
            updated_after=updated_after,
        )
        profiles = self._fetch_profiles(
            database,
            scheme_id,
            added_after=added_after,
            updated_after=updated_after,
        )

        return SchemeData(
            metadata=metadata,
            loci=loci,
            profiles=profiles,
            cutoff_date=cutoff_date,
        )
