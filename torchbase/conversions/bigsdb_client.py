"""BIGSdb REST API client for scheme retrieval.

This module provides a client for accessing BIGSdb REST endpoints to retrieve
scheme metadata, loci information, allele sequences, and profile tables.
It handles automatic pagination, temporal filtering, and captures provenance
timestamps for all responses.
"""

from dataclasses import dataclass
from datetime import date, datetime, timezone
from typing import List, Optional, Dict, Any, Tuple
import csv
import io
import re
import requests


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

    def __init__(self, base_url: str, timeout: int = 30):
        """Initialize BIGSdb client.

        Args:
            base_url: Base URL of BIGSdb REST API
            timeout: Request timeout in seconds
        """
        self.base_url = base_url
        self.timeout = timeout

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

        try:
            response = requests.get(url, params=params, timeout=self.timeout)
        except (
            requests.ConnectionError,
            requests.Timeout,
            requests.RequestException,
            ConnectionError,
            TimeoutError,
        ) as e:
            raise BIGSdbNetworkError(f"Network error: {str(e)}") from e

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

        if not data.get("records") or len(data["records"]) == 0:
            msg = f"Scheme {scheme_id} not found in database {database}"
            raise BIGSdbError(msg)

        record = data["records"][0]
        last_updated = None
        if record.get("last_updated"):
            last_updated = self._parse_datetime(record["last_updated"])

        return SchemeMetadata(
            scheme_id=int(record.get("scheme_id", scheme_id)),
            name=record.get("description", f"Scheme {scheme_id}"),
            description=record.get("description"),
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

        while True:
            params["page"] = page
            endpoint = f"db/{database}/schemes/{scheme_id}/loci"
            data = self._make_request("GET", endpoint, params=params)

            records = data.get("records", [])
            for record in records:
                last_updated = None
                if record.get("last_updated"):
                    last_updated = self._parse_datetime(record["last_updated"])

                locus = LocusData(
                    locus_id=record.get("id", ""),
                    locus_name=record.get("locus", ""),
                    alleles_count=int(record.get("alleles", 0)),
                    last_updated=last_updated,
                )
                loci.append(locus)

            paging = data.get("paging", {})
            total_pages = paging.get("pages", 1)
            if page >= total_pages:
                break

            page += 1

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

        # Capture current timestamp for provenance
        last_updated = datetime.now(timezone.utc)

        return ProfileTable(
            profiles=profiles,
            row_count=len(profiles),
            last_updated=last_updated,
        )

    def _fetch_alleles_fasta(
        self, database: str, locus_id: str, cutoff_date: Optional[date] = None
    ) -> str:
        """Fetch allele sequences for a locus as FASTA text.

        When cutoff_date is given, fetches via the JSON alleles endpoint and
        filters by date_entered to stay within freely-licensed data (pre-2025).
        Without a cutoff, hits the alleles_fasta endpoint directly.

        Args:
            database: Database identifier (e.g., "pubmlst")
            locus_id: Locus identifier string
            cutoff_date: Exclude alleles entered after this date

        Returns:
            FASTA text containing allele sequences

        Raises:
            BIGSdbError: If alleles cannot be fetched
        """
        if cutoff_date is not None:
            return self._fetch_alleles_fasta_filtered(database, locus_id, cutoff_date)
        endpoint = f"db/{database}/loci/{locus_id}/alleles_fasta"
        return self._make_request("GET", endpoint, expect_json=False)

    def _fetch_alleles_fasta_filtered(
        self, database: str, locus_id: str, cutoff_date: date
    ) -> str:
        """Fetch alleles via JSON endpoint, filter by cutoff_date, return as FASTA.

        PubMLST changed their licensing on 2025-01-01: alleles submitted from
        that date onward require authentication and carry non-commercial terms.
        This method fetches individual allele records (which include date_entered)
        and excludes any entered after the cutoff, producing a FASTA of only the
        freely-distributable sequences.

        Args:
            database: Database identifier
            locus_id: Locus identifier
            cutoff_date: Exclude alleles with date_entered after this date

        Returns:
            FASTA text of filtered alleles

        Raises:
            BIGSdbError: If alleles cannot be fetched
        """
        lines = []
        page = 1

        while True:
            endpoint = f"db/{database}/loci/{locus_id}/alleles"
            data = self._make_request("GET", endpoint, params={"page": page})

            records = data.get("alleles", data.get("records", []))
            for record in records:
                date_str = record.get("date_entered", "")
                if date_str:
                    try:
                        entry_date = datetime.strptime(date_str[:10], "%Y-%m-%d").date()
                        if entry_date > cutoff_date:
                            continue
                    except ValueError:
                        pass

                allele_id = record.get("allele_id", "")
                sequence = record.get("sequence", "")
                if allele_id and sequence:
                    lines.append(f">{locus_id}_{allele_id}")
                    lines.append(sequence)

            paging = data.get("paging", {})
            total_pages = paging.get("pages", 1)
            if page >= total_pages:
                break
            page += 1

        return "\n".join(lines) + ("\n" if lines else "")

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
        data = self._make_request("GET", "databases")
        raw = data.get("databases", data.get("records", []))

        databases = []
        for entry in raw:
            # entry may be a dict or a URL string
            if isinstance(entry, str):
                name = entry.rstrip("/").split("/")[-1]
                description = name
            else:
                name = entry.get("name", entry.get("href", "").rstrip("/").split("/")[-1])
                description = entry.get("description", name)

            if name.endswith("_seqdef"):
                databases.append(DatabaseInfo(name=name, description=description))

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
