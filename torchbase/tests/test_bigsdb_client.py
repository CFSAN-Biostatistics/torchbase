"""Tests for BIGSdb REST client."""

import pytest
from datetime import datetime
from unittest.mock import Mock, patch
import json

from torchbase.conversions.bigsdb_client import (
    BIGSdbClient,
    SchemeMetadata,
    LocusData,
    ProfileTable,
    SchemeData,
    BIGSdbError,
    BIGSdbNetworkError,
    BIGSdbValidationError,
)


class TestSchemeMetadataDataclass:
    """Test SchemeMetadata dataclass"""

    def test_create_scheme_metadata(self):
        metadata = SchemeMetadata(
            scheme_id=1,
            name="MLST",
            description="Multi-locus sequence typing",
            last_updated=datetime(2023, 1, 15, 10, 30, 0),
        )
        assert metadata.scheme_id == 1
        assert metadata.name == "MLST"
        assert metadata.description == "Multi-locus sequence typing"
        assert metadata.last_updated == datetime(2023, 1, 15, 10, 30, 0)


class TestLocusDataDataclass:
    """Test LocusData dataclass"""

    def test_create_locus_data(self):
        locus = LocusData(
            locus_id="adk",
            locus_name="Adenylate kinase",
            alleles_count=100,
            last_updated=datetime(2023, 1, 10, 9, 0, 0),
        )
        assert locus.locus_id == "adk"
        assert locus.locus_name == "Adenylate kinase"
        assert locus.alleles_count == 100
        assert locus.last_updated == datetime(2023, 1, 10, 9, 0, 0)


class TestProfileTableDataclass:
    """Test ProfileTable dataclass"""

    def test_create_profile_table(self):
        table = ProfileTable(
            profiles=[
                {"ST": "1", "adk": "1", "fumC": "1", "gyrB": "1"},
                {"ST": "2", "adk": "2", "fumC": "2", "gyrB": "2"},
            ],
            row_count=2,
            last_updated=datetime(2023, 1, 12, 8, 0, 0),
        )
        assert len(table.profiles) == 2
        assert table.profiles[0]["ST"] == "1"
        assert table.row_count == 2
        assert table.last_updated == datetime(2023, 1, 12, 8, 0, 0)


class TestSchemeDataDataclass:
    """Test SchemeData dataclass"""

    def test_create_scheme_data(self):
        metadata = SchemeMetadata(
            scheme_id=1,
            name="MLST",
            description="Multi-locus sequence typing",
            last_updated=datetime(2023, 1, 15, 10, 30, 0),
        )
        loci = [
            LocusData(
                locus_id="adk",
                locus_name="Adenylate kinase",
                alleles_count=100,
                last_updated=datetime(2023, 1, 10, 9, 0, 0),
            )
        ]
        scheme_data = SchemeData(
            metadata=metadata,
            loci=loci,
        )
        assert scheme_data.metadata == metadata
        assert len(scheme_data.loci) == 1
        assert scheme_data.loci[0].locus_id == "adk"


class TestBIGSdbClientInit:
    """Test BIGSdbClient initialization"""

    def test_init_with_base_url(self):
        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        assert client.base_url == "http://pubmlst.org/api"

    def test_init_with_default_timeout(self):
        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        assert client.timeout == 30

    def test_init_with_custom_timeout(self):
        client = BIGSdbClient(base_url="http://pubmlst.org/api", timeout=60)
        assert client.timeout == 60


class TestBIGSdbClientFetchSchemeMetadata:
    """Test fetching scheme metadata"""

    @patch("torchbase.conversions.bigsdb_client.requests.get")
    def test_fetch_scheme_metadata_success(self, mock_get):
        """Test successful fetch of scheme metadata"""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "records": [
                {
                    "scheme_id": "1",
                    "description": "MLST",
                    "last_updated": "2023-01-15T10:30:00",
                }
            ]
        }
        mock_get.return_value = mock_response

        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        metadata = client._fetch_scheme_metadata("pubmlst", 1)

        assert metadata.scheme_id == 1
        assert metadata.name == "MLST"
        assert metadata.last_updated == datetime(2023, 1, 15, 10, 30, 0)
        mock_get.assert_called_once()

    @patch("torchbase.conversions.bigsdb_client.requests.get")
    def test_fetch_scheme_metadata_http_error(self, mock_get):
        """Test HTTP error handling"""
        mock_response = Mock()
        mock_response.status_code = 500
        mock_get.return_value = mock_response

        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        with pytest.raises(BIGSdbError):
            client._fetch_scheme_metadata("pubmlst", 1)

    @patch("torchbase.conversions.bigsdb_client.requests.get")
    def test_fetch_scheme_metadata_not_found(self, mock_get):
        """Test 404 for invalid scheme"""
        mock_response = Mock()
        mock_response.status_code = 404
        mock_get.return_value = mock_response

        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        with pytest.raises(BIGSdbError):
            client._fetch_scheme_metadata("pubmlst", 999)


class TestBIGSdbClientFetchLoci:
    """Test fetching loci data"""

    @patch("torchbase.conversions.bigsdb_client.requests.get")
    def test_fetch_loci_success_single_page(self, mock_get):
        """Test successful fetch of loci without pagination"""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "records": [
                {
                    "id": "adk",
                    "locus": "Adenylate kinase",
                    "alleles": "100",
                    "last_updated": "2023-01-10T09:00:00",
                }
            ],
            "paging": {"pages": 1, "return_limit": 1},
        }
        mock_get.return_value = mock_response

        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        loci = client._fetch_loci("pubmlst", 1)

        assert len(loci) == 1
        assert loci[0].locus_id == "adk"
        assert loci[0].alleles_count == 100

    @patch("torchbase.conversions.bigsdb_client.requests.get")
    def test_fetch_loci_with_pagination(self, mock_get):
        """Test automatic pagination handling"""
        # First page
        first_response = Mock()
        first_response.status_code = 200
        first_response.json.return_value = {
            "records": [
                {
                    "id": "adk",
                    "locus": "Adenylate kinase",
                    "alleles": "100",
                    "last_updated": "2023-01-10T09:00:00",
                }
            ],
            "paging": {"pages": 2, "return_limit": 1},
        }

        # Second page
        second_response = Mock()
        second_response.status_code = 200
        second_response.json.return_value = {
            "records": [
                {
                    "id": "fumC",
                    "locus": "Fumarate hydratase",
                    "alleles": "95",
                    "last_updated": "2023-01-11T09:00:00",
                }
            ],
            "paging": {"pages": 2, "return_limit": 1},
        }

        mock_get.side_effect = [first_response, second_response]

        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        loci = client._fetch_loci("pubmlst", 1)

        assert len(loci) == 2
        assert loci[0].locus_id == "adk"
        assert loci[1].locus_id == "fumC"
        assert mock_get.call_count == 2

    @patch("torchbase.conversions.bigsdb_client.requests.get")
    def test_fetch_loci_with_added_after_filter(self, mock_get):
        """Test temporal filtering with added_after"""  # noqa: E501
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "records": [
                {
                    "id": "adk",
                    "locus": "Adenylate kinase",
                    "alleles": "100",
                    "last_updated": "2023-01-10T09:00:00",
                }
            ],
            "paging": {"pages": 1, "return_limit": 1},
        }
        mock_get.return_value = mock_response

        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        loci = client._fetch_loci(
            "pubmlst", 1, added_after=datetime(2023, 1, 1, 0, 0, 0)
        )

        assert len(loci) == 1
        # Verify that the filter was passed in the URL
        call_args = mock_get.call_args
        assert "added_after" in call_args[1].get("params", {})


class TestBIGSdbClientFetchProfiles:
    """Test fetching profile data"""

    @patch("torchbase.conversions.bigsdb_client.requests.get")
    def test_fetch_profiles_success(self, mock_get):
        """Test successful fetch of profiles as CSV"""
        csv_content = "ST\tadk\tfumC\tgyrB\n1\t1\t1\t1\n2\t2\t2\t2\n"
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = csv_content
        mock_get.return_value = mock_response

        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        profiles = client._fetch_profiles("pubmlst", 1)

        assert isinstance(profiles, ProfileTable)
        assert len(profiles.profiles) == 2
        assert profiles.profiles[0]["ST"] == "1"
        assert profiles.profiles[1]["ST"] == "2"

    @patch("torchbase.conversions.bigsdb_client.requests.get")
    def test_fetch_profiles_with_updated_after_filter(self, mock_get):
        """Test temporal filtering with updated_after"""  # noqa: E501
        csv_content = "ST\tadk\tfumC\tgyrB\n1\t1\t1\t1\n"
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = csv_content
        mock_get.return_value = mock_response

        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        profiles = client._fetch_profiles(
            "pubmlst", 1, updated_after=datetime(2023, 1, 1, 0, 0, 0)
        )

        assert len(profiles.profiles) == 1
        # Verify that the filter was passed in the URL
        call_args = mock_get.call_args
        assert "updated_after" in call_args[1].get("params", {})


class TestBIGSdbClientFetchScheme:
    """Test fetch_scheme - main public API"""

    @patch("torchbase.conversions.bigsdb_client.requests.get")
    def test_fetch_scheme_success(self, mock_get):
        """Test successful scheme fetch with all components"""
        # Mock scheme metadata response
        scheme_response = Mock()
        scheme_response.status_code = 200
        scheme_response.json.return_value = {
            "records": [
                {
                    "scheme_id": "1",
                    "description": "MLST",
                    "last_updated": "2023-01-15T10:30:00",
                }
            ]
        }

        # Mock loci response
        loci_response = Mock()
        loci_response.status_code = 200
        loci_response.json.return_value = {
            "records": [
                {
                    "id": "adk",
                    "locus": "Adenylate kinase",
                    "alleles": "100",
                    "last_updated": "2023-01-10T09:00:00",
                }
            ],
            "paging": {"pages": 1, "return_limit": 1},
        }

        # Mock profiles response
        csv_content = "ST\tadk\n1\t1\n2\t2\n"
        profiles_response = Mock()
        profiles_response.status_code = 200
        profiles_response.text = csv_content

        mock_get.side_effect = [
            scheme_response,
            loci_response,
            profiles_response,
        ]

        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        scheme_data = client.fetch_scheme("pubmlst", 1)

        assert isinstance(scheme_data, SchemeData)
        assert scheme_data.metadata.scheme_id == 1
        assert scheme_data.metadata.name == "MLST"
        assert len(scheme_data.loci) == 1
        assert len(scheme_data.profiles.profiles) == 2

    @patch("torchbase.conversions.bigsdb_client.requests.get")
    def test_fetch_scheme_with_filters(self, mock_get):
        """Test fetch_scheme with temporal filters"""
        scheme_response = Mock()
        scheme_response.status_code = 200
        scheme_response.json.return_value = {
            "records": [
                {
                    "scheme_id": "1",
                    "description": "MLST",
                    "last_updated": "2023-01-15T10:30:00",
                }
            ]
        }

        loci_response = Mock()
        loci_response.status_code = 200
        loci_response.json.return_value = {
            "records": [],
            "paging": {"pages": 1, "return_limit": 0},
        }

        csv_content = "ST\tadk\n"
        profiles_response = Mock()
        profiles_response.status_code = 200
        profiles_response.text = csv_content

        mock_get.side_effect = [
            scheme_response,
            loci_response,
            profiles_response,
        ]

        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        scheme_data = client.fetch_scheme(
            "pubmlst",
            1,
            added_after=datetime(2023, 1, 1, 0, 0, 0),
            updated_after=datetime(2023, 1, 1, 0, 0, 0),
        )

        assert isinstance(scheme_data, SchemeData)


class TestBIGSdbErrorHandling:
    """Test error handling and exceptions"""

    def test_bigsdb_error_exception(self):
        """Test BIGSdbError exception"""
        error = BIGSdbError("Test error")
        assert str(error) == "Test error"

    def test_bigsdb_network_error_exception(self):
        """Test BIGSdbNetworkError exception"""
        error = BIGSdbNetworkError("Connection failed")
        assert str(error) == "Connection failed"

    def test_bigsdb_validation_error_exception(self):
        """Test BIGSdbValidationError exception"""
        error = BIGSdbValidationError("Invalid data")
        assert str(error) == "Invalid data"

    @patch("torchbase.conversions.bigsdb_client.requests.get")
    def test_network_error_handling(self, mock_get):
        """Test handling of network errors"""  # noqa: E501
        mock_get.side_effect = ConnectionError("Network timeout")

        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        with pytest.raises(BIGSdbNetworkError):
            client._fetch_scheme_metadata("pubmlst", 1)  # noqa: E501

    @patch("torchbase.conversions.bigsdb_client.requests.get")
    def test_json_decode_error_handling(self, mock_get):
        """Test handling of invalid JSON responses"""
        mock_response = Mock()
        mock_response.status_code = 200
        err = json.JSONDecodeError("Invalid JSON", "", 0)
        mock_response.json.side_effect = err
        mock_get.return_value = mock_response

        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        with pytest.raises(BIGSdbValidationError):
            client._fetch_scheme_metadata("pubmlst", 1)


class TestBIGSdbTimestampCapture:
    """Test timestamp capture for provenance"""

    @patch("torchbase.conversions.bigsdb_client.requests.get")
    def test_scheme_metadata_captures_timestamp(self, mock_get):
        """Test that scheme metadata captures last_updated"""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "records": [
                {
                    "scheme_id": "1",
                    "description": "MLST",
                    "last_updated": "2023-01-15T10:30:00",
                }
            ]
        }
        mock_get.return_value = mock_response

        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        metadata = client._fetch_scheme_metadata("pubmlst", 1)

        assert metadata.last_updated is not None
        assert isinstance(metadata.last_updated, datetime)

    @patch("torchbase.conversions.bigsdb_client.requests.get")
    def test_locus_data_captures_timestamp(self, mock_get):
        """Test that locus data captures last_updated"""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "records": [
                {
                    "id": "adk",
                    "locus": "Adenylate kinase",
                    "alleles": "100",
                    "last_updated": "2023-01-10T09:00:00",
                }
            ],
            "paging": {"pages": 1, "return_limit": 1},
        }
        mock_get.return_value = mock_response

        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        loci = client._fetch_loci("pubmlst", 1)

        assert loci[0].last_updated is not None
        assert isinstance(loci[0].last_updated, datetime)

    @patch("torchbase.conversions.bigsdb_client.requests.get")
    def test_profile_table_captures_timestamp(self, mock_get):
        """Test that profile table captures last_updated"""
        csv_content = "ST\tadk\n1\t1\n"
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = csv_content
        mock_get.return_value = mock_response

        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        profiles = client._fetch_profiles("pubmlst", 1)

        assert profiles.last_updated is not None
        assert isinstance(profiles.last_updated, datetime)
