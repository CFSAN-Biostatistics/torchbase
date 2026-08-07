"""Tests for BIGSdb REST client."""

import pytest
from datetime import datetime
from pathlib import Path
from unittest.mock import Mock, MagicMock, patch, call
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


def _make_client(base_url="http://pubmlst.org/api", mock_session=None, **kwargs):
    """Create a BIGSdbClient with an injected mock session for testing."""
    client = BIGSdbClient(base_url=base_url, **kwargs)
    if mock_session is not None:
        client._session = mock_session
    return client


def _mock_response(status_code=200, json_data=None, text=None):
    resp = Mock()
    resp.status_code = status_code
    resp.reason = "OK" if status_code < 400 else "Error"
    if json_data is not None:
        resp.json.return_value = json_data
    if text is not None:
        resp.text = text
    resp.headers = {}
    resp.iter_content = Mock(return_value=iter([]))
    return resp


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

    def test_init_default_verify_true(self):
        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        assert client.verify is True

    def test_init_verify_false(self):
        client = BIGSdbClient(base_url="http://pubmlst.org/api", verify=False)
        assert client.verify is False

    def test_init_verify_ca_bundle_path(self):
        client = BIGSdbClient(base_url="http://pubmlst.org/api", verify="/path/to/bundle.pem")
        assert client.verify == "/path/to/bundle.pem"

    def test_init_creates_session(self):
        import requests
        client = BIGSdbClient(base_url="http://pubmlst.org/api")
        assert isinstance(client._session, requests.Session)


class TestBIGSdbClientFetchSchemeMetadata:
    """Test fetching scheme metadata"""

    def test_fetch_scheme_metadata_success(self):
        mock_session = Mock()
        mock_session.get.return_value = _mock_response(json_data={
            "id": 1,
            "description": "MLST",
            "last_updated": "2023-01-15",
        })

        client = _make_client(mock_session=mock_session)
        metadata = client._fetch_scheme_metadata("pubmlst", 1)

        assert metadata.scheme_id == 1
        assert metadata.name == "MLST"
        assert metadata.last_updated is not None
        mock_session.get.assert_called_once()

    def test_fetch_scheme_metadata_http_error(self):
        mock_session = Mock()
        mock_session.get.return_value = _mock_response(status_code=500)

        client = _make_client(mock_session=mock_session)
        with pytest.raises(BIGSdbError):
            client._fetch_scheme_metadata("pubmlst", 1)

    def test_fetch_scheme_metadata_not_found(self):
        mock_session = Mock()
        mock_session.get.return_value = _mock_response(status_code=404)

        client = _make_client(mock_session=mock_session)
        with pytest.raises(BIGSdbError):
            client._fetch_scheme_metadata("pubmlst", 999)

    def test_verify_propagated_to_request(self):
        mock_session = Mock()
        mock_session.get.return_value = _mock_response(json_data={
            "id": 1, "description": "MLST",
        })

        client = _make_client(mock_session=mock_session, verify=False)
        client._fetch_scheme_metadata("pubmlst", 1)

        call_kwargs = mock_session.get.call_args[1]
        assert call_kwargs.get("verify") is False

    def test_ca_bundle_propagated_to_request(self):
        mock_session = Mock()
        mock_session.get.return_value = _mock_response(json_data={
            "id": 1, "description": "MLST",
        })

        client = _make_client(mock_session=mock_session, verify="/path/to/bundle.pem")
        client._fetch_scheme_metadata("pubmlst", 1)

        call_kwargs = mock_session.get.call_args[1]
        assert call_kwargs.get("verify") == "/path/to/bundle.pem"


class TestBIGSdbClientFetchLoci:
    """Test fetching loci data"""

    def test_fetch_loci_success_single_page(self):
        mock_session = Mock()
        mock_session.get.return_value = _mock_response(json_data={
            "loci": ["http://pubmlst.org/api/db/pubmlst/loci/adk"],
            "records": 1,
        })

        client = _make_client(mock_session=mock_session)
        loci = client._fetch_loci("pubmlst", 1)

        assert len(loci) == 1
        assert loci[0].locus_id == "adk"

    def test_fetch_loci_with_pagination(self):
        mock_session = Mock()
        mock_session.get.return_value = _mock_response(json_data={
            "loci": [
                "http://pubmlst.org/api/db/pubmlst/loci/adk",
                "http://pubmlst.org/api/db/pubmlst/loci/fumC",
            ],
            "records": 2,
        })

        client = _make_client(mock_session=mock_session)
        loci = client._fetch_loci("pubmlst", 1)

        assert len(loci) == 2
        assert loci[0].locus_id == "adk"
        assert loci[1].locus_id == "fumC"

    def test_fetch_loci_with_added_after_filter(self):
        mock_session = Mock()
        mock_session.get.return_value = _mock_response(json_data={
            "loci": ["http://pubmlst.org/api/db/pubmlst/loci/adk"],
            "records": 1,
        })

        client = _make_client(mock_session=mock_session)
        loci = client._fetch_loci(
            "pubmlst", 1, added_after=datetime(2023, 1, 1, 0, 0, 0)
        )

        assert len(loci) == 1
        call_kwargs = mock_session.get.call_args[1]
        assert "added_after" in call_kwargs.get("params", {})


class TestBIGSdbClientFetchProfiles:
    """Test fetching profile data"""

    def test_fetch_profiles_success(self):
        csv_content = "ST\tadk\tfumC\tgyrB\n1\t1\t1\t1\n2\t2\t2\t2\n"
        mock_session = Mock()
        mock_session.get.return_value = _mock_response(status_code=200, text=csv_content)

        client = _make_client(mock_session=mock_session)
        profiles = client._fetch_profiles("pubmlst", 1)

        assert isinstance(profiles, ProfileTable)
        assert len(profiles.profiles) == 2
        assert profiles.profiles[0]["ST"] == "1"
        assert profiles.profiles[1]["ST"] == "2"

    def test_fetch_profiles_with_updated_after_filter(self):
        csv_content = "ST\tadk\tfumC\tgyrB\n1\t1\t1\t1\n"
        mock_session = Mock()
        mock_session.get.return_value = _mock_response(status_code=200, text=csv_content)

        client = _make_client(mock_session=mock_session)
        profiles = client._fetch_profiles(
            "pubmlst", 1, updated_after=datetime(2023, 1, 1, 0, 0, 0)
        )

        assert len(profiles.profiles) == 1
        call_kwargs = mock_session.get.call_args[1]
        assert "updated_after" in call_kwargs.get("params", {})


class TestBIGSdbClientFetchAllelesStreaming:
    """Test _fetch_alleles_fasta streaming-to-file behaviour"""

    def test_fetch_alleles_writes_to_dest(self, tmp_path):
        fasta_bytes = b">allele_1\nACGT\n>allele_2\nTGCA\n"
        mock_session = Mock()
        resp = Mock()
        resp.status_code = 200
        resp.headers = {"content-length": str(len(fasta_bytes))}
        resp.iter_content.return_value = iter([fasta_bytes])
        mock_session.get.return_value = resp

        client = _make_client(mock_session=mock_session)
        dest = tmp_path / "adk.fasta"
        result = client._fetch_alleles_fasta("pubmlst_neisseria_seqdef", "adk", dest)

        assert result == dest
        assert dest.exists()
        assert dest.read_bytes() == fasta_bytes

    def test_fetch_alleles_raises_on_http_error(self, tmp_path):
        mock_session = Mock()
        resp = Mock()
        resp.status_code = 404
        resp.reason = "Not Found"
        mock_session.get.return_value = resp

        client = _make_client(mock_session=mock_session)
        with pytest.raises(BIGSdbError):
            client._fetch_alleles_fasta("pubmlst_neisseria_seqdef", "adk", tmp_path / "adk.fasta")

    def test_fetch_alleles_uses_stream_true(self, tmp_path):
        mock_session = Mock()
        resp = Mock()
        resp.status_code = 200
        resp.headers = {}
        resp.iter_content.return_value = iter([b">a\nACGT\n"])
        mock_session.get.return_value = resp

        client = _make_client(mock_session=mock_session)
        client._fetch_alleles_fasta("pubmlst_neisseria_seqdef", "adk", tmp_path / "adk.fasta")

        call_kwargs = mock_session.get.call_args[1]
        assert call_kwargs.get("stream") is True

    def test_fetch_alleles_verify_propagated(self, tmp_path):
        mock_session = Mock()
        resp = Mock()
        resp.status_code = 200
        resp.headers = {}
        resp.iter_content.return_value = iter([b">a\nACGT\n"])
        mock_session.get.return_value = resp

        client = _make_client(mock_session=mock_session, verify=False)
        client._fetch_alleles_fasta("pubmlst_neisseria_seqdef", "adk", tmp_path / "adk.fasta")

        call_kwargs = mock_session.get.call_args[1]
        assert call_kwargs.get("verify") is False

    def test_fetch_alleles_creates_parent_dirs(self, tmp_path):
        mock_session = Mock()
        resp = Mock()
        resp.status_code = 200
        resp.headers = {}
        resp.iter_content.return_value = iter([b">a\nACGT\n"])
        mock_session.get.return_value = resp

        client = _make_client(mock_session=mock_session)
        dest = tmp_path / "deep" / "nested" / "adk.fasta"
        client._fetch_alleles_fasta("pubmlst_neisseria_seqdef", "adk", dest)
        assert dest.exists()

    def test_fetch_alleles_network_error_raises(self, tmp_path):
        import requests as _requests
        mock_session = Mock()
        mock_session.get.side_effect = _requests.ConnectionError("no route")

        client = _make_client(mock_session=mock_session)
        with pytest.raises(BIGSdbNetworkError):
            client._fetch_alleles_fasta("pubmlst_neisseria_seqdef", "adk", tmp_path / "adk.fasta")


class TestBIGSdbClientFetchScheme:
    """Test fetch_scheme - main public API"""

    def test_fetch_scheme_success(self):
        scheme_response = _mock_response(json_data={
            "id": 1, "description": "MLST", "last_updated": "2023-01-15",
        })
        loci_response = _mock_response(json_data={
            "loci": ["http://pubmlst.org/api/db/pubmlst/loci/adk"],
            "records": 1,
        })
        csv_content = "ST\tadk\n1\t1\n2\t2\n"
        profiles_response = _mock_response(status_code=200, text=csv_content)

        mock_session = Mock()
        mock_session.get.side_effect = [scheme_response, loci_response, profiles_response]

        client = _make_client(mock_session=mock_session)
        scheme_data = client.fetch_scheme("pubmlst", 1)

        assert isinstance(scheme_data, SchemeData)
        assert scheme_data.metadata.scheme_id == 1
        assert scheme_data.metadata.name == "MLST"
        assert len(scheme_data.loci) == 1
        assert len(scheme_data.profiles.profiles) == 2

    def test_fetch_scheme_with_filters(self):
        mock_session = Mock()
        mock_session.get.side_effect = [
            _mock_response(json_data={"id": 1, "description": "MLST", "last_updated": "2023-01-15"}),
            _mock_response(json_data={"loci": [], "records": 0}),
            _mock_response(status_code=200, text="ST\tadk\n"),
        ]

        client = _make_client(mock_session=mock_session)
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
        error = BIGSdbError("Test error")
        assert str(error) == "Test error"

    def test_bigsdb_network_error_exception(self):
        error = BIGSdbNetworkError("Connection failed")
        assert str(error) == "Connection failed"

    def test_bigsdb_validation_error_exception(self):
        error = BIGSdbValidationError("Invalid data")
        assert str(error) == "Invalid data"

    def test_network_error_handling(self):
        import requests as _requests
        mock_session = Mock()
        mock_session.get.side_effect = ConnectionError("Network timeout")

        client = _make_client(mock_session=mock_session)
        with pytest.raises(BIGSdbNetworkError):
            client._fetch_scheme_metadata("pubmlst", 1)

    def test_json_decode_error_handling(self):
        mock_session = Mock()
        resp = Mock()
        resp.status_code = 200
        resp.json.side_effect = json.JSONDecodeError("Invalid JSON", "", 0)
        mock_session.get.return_value = resp

        client = _make_client(mock_session=mock_session)
        with pytest.raises(BIGSdbValidationError):
            client._fetch_scheme_metadata("pubmlst", 1)


class TestBIGSdbTimestampCapture:
    """Test timestamp capture for provenance"""

    def test_scheme_metadata_captures_timestamp(self):
        mock_session = Mock()
        mock_session.get.return_value = _mock_response(json_data={
            "id": 1, "description": "MLST", "last_updated": "2023-01-15",
        })

        client = _make_client(mock_session=mock_session)
        metadata = client._fetch_scheme_metadata("pubmlst", 1)

        assert metadata.last_updated is not None
        assert isinstance(metadata.last_updated, datetime)

    def test_locus_data_is_returned(self):
        mock_session = Mock()
        mock_session.get.return_value = _mock_response(json_data={
            "loci": ["http://pubmlst.org/api/db/pubmlst/loci/adk"],
            "records": 1,
        })

        client = _make_client(mock_session=mock_session)
        loci = client._fetch_loci("pubmlst", 1)

        assert len(loci) == 1
        assert loci[0].locus_id == "adk"

    def test_profile_table_captures_timestamp(self):
        mock_session = Mock()
        mock_session.get.return_value = _mock_response(status_code=200, text="ST\tadk\n1\t1\n")

        client = _make_client(mock_session=mock_session)
        profiles = client._fetch_profiles("pubmlst", 1)

        assert profiles.last_updated is not None
        assert isinstance(profiles.last_updated, datetime)
