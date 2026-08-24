"""Tests for torchbase.typing_run's calling-mode dispatch (docs/adr/0003).

A torch with `calling_mode == "presence_absence"` always types at
`sensitive` (alignment-only) and reduces alignments through
`allele_calls.calls_from_presence` instead of `calls_from_alignment`, with
`id_column` threaded into profile matching. Compute (minimap2/WDL) is
stubbed throughout: only the dispatch decision is under test.
"""

from unittest.mock import MagicMock, patch

import pytest

from torchbase import typing_run


def _torch(tmp_path, alleles, profiles_tsv, calling_mode="identity", id_column=None):
    torch = MagicMock()
    allele_fasta = tmp_path / "alleles.fasta"
    allele_fasta.write_text(alleles)
    profiles = tmp_path / "profiles.tsv"
    profiles.write_text(profiles_tsv)
    torch.get_unified_files.return_value = (allele_fasta, profiles)
    torch.path = tmp_path / "ns" / "name" / "1.0.0.torch"
    torch.calling_mode = calling_mode
    torch.id_column = id_column
    return torch


def _query(tmp_path):
    query = tmp_path / "query.fasta"
    query.write_text(">contig_1\nACGT\n")
    return query


class TestPresenceAbsenceDispatch:
    def test_presence_absence_torch_never_screens(self, tmp_path):
        """MinHash screening is skipped entirely -- alignment always runs."""
        torch = _torch(
            tmp_path, ">cadA\nACGT\n", "Serotype\tcadA\nX\tpresent\n",
            calling_mode="presence_absence", id_column="Serotype",
        )
        screen_called = []
        with patch.object(typing_run, "_screen", side_effect=lambda *a, **k: screen_called.append(1)), \
             patch.object(typing_run, "_align_presence", return_value={"cadA": {"allele_id": "present"}}):
            typing_run.type_allelic(torch, str(_query(tmp_path)), strategy="fast")
        assert screen_called == []

    def test_presence_absence_torch_ignores_requested_strategy(self, tmp_path):
        """Even an explicit fast/balanced strategy is overridden to sensitive."""
        torch = _torch(
            tmp_path, ">cadA\nACGT\n", "Serotype\tcadA\nX\tpresent\n",
            calling_mode="presence_absence", id_column="Serotype",
        )
        with patch.object(typing_run, "_align_presence", return_value={"cadA": {"allele_id": "present"}}) as mock_align:
            result = typing_run.type_allelic(torch, str(_query(tmp_path)), strategy="balanced")
        mock_align.assert_called_once()
        assert result["method"]["strategy"] == "sensitive"

    def test_presence_absence_calls_reach_the_matcher_with_id_column(self, tmp_path):
        torch = _torch(
            tmp_path, ">cadA\nACGT\n",
            "Serotype\tcadA\nSonnei\tpresent\n",
            calling_mode="presence_absence", id_column="Serotype",
        )
        with patch.object(typing_run, "_align_presence", return_value={"cadA": {"allele_id": "present"}}):
            result = typing_run.type_allelic(torch, str(_query(tmp_path)), strategy="balanced")
        assert result["profile_id"] == "Sonnei"
        assert result["status"] == "known_profile"

    def test_identity_torch_is_unaffected(self, tmp_path):
        """The default calling mode still uses the screen/align split as before."""
        torch = _torch(
            tmp_path, ">adk_1\nACGT\n", "ST\tadk\n1\t1\n",
            calling_mode="identity", id_column=None,
        )
        with patch.object(
            typing_run, "_screen",
            return_value={"adk": {"allele_id": "1", "confidence": True}},
        ) as mock_screen, patch.object(typing_run, "_align_presence") as mock_align_presence:
            result = typing_run.type_allelic(torch, str(_query(tmp_path)), strategy="fast")
        mock_screen.assert_called_once()
        mock_align_presence.assert_not_called()
        assert result["profile_id"] == "1"
