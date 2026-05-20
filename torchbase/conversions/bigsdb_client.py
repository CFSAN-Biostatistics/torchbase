"""BIGSdb REST API client for scheme retrieval.

This module provides a client for accessing BIGSdb REST endpoints to retrieve
scheme metadata, loci information, allele sequences, and profile tables.
It handles automatic pagination, temporal filtering, and captures provenance
timestamps for all responses.
"""

from dataclasses import dataclass
from datetime import datetime, timezone
from typing import List, Optional, Dict, Any
import csv
import io
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
    """

    metadata: SchemeMetadata
    loci: List[LocusData]
    profiles: Optional[ProfileTable] = None


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

    def fetch_scheme(
        self,
        database: str,
        scheme_id: int,
        added_after: Optional[datetime] = None,
        updated_after: Optional[datetime] = None,
    ) -> SchemeData:
        """Fetch complete scheme data.

        Retrieves scheme metadata, loci, and profiles in one call.
        Automatically handles pagination for loci and applies
        temporal filters.

        Args:
            database: Database identifier (e.g., "pubmlst")
            scheme_id: Numeric scheme identifier
            added_after: Filter for records added after this datetime
            updated_after: Filter for records updated after this
                datetime

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
        )
