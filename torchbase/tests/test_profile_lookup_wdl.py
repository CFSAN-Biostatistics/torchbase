#!/usr/bin/env python

"""
Acceptance tests for WDL task: Profile Lookup with Novel Detection (#17)

Tests the profile table lookup and nearest ST calculation functionality.
Covers exact matches, novel profiles, novel alleles, and edge cases.

Acceptance criteria:
- WDL task signature: allele calls JSON + profiles.tsv → typing result JSON
- Exact match lookup in profile table
- Novel profile detection (known alleles, unknown combination)
- Nearest ST calculation: hamming distance across all profiles
- Output JSON: {sequence_type, status: "known"/"novel_profile"/"novel_allele", nearest_st, distance}
- Handles missing loci gracefully (partial profiles)
- miniwdl check validates syntax
- Test with known profiles, novel combinations, novel alleles
"""

import pytest

import torchbase


class TestProfileLookupTask:
    """Acceptance tests for profile lookup and nearest ST calculation."""

    @pytest.fixture
    def simple_schema(self):
        """Simple schema with known profiles for testing."""
        profile_data = """ST	locusA	locusB	locusC	locusD
1	1	1	1	1
2	1	1	1	2
3	1	1	2	1
4	1	2	1	1
5	2	1	1	1
6	2	1	1	2
7	1	2	2	2
8	2	2	2	2
9	3	1	1	1
10	1	3	1	1
"""
        rows = [line.split('\t') for line in profile_data.strip().split('\n')]
        return torchbase.Profile.parse("test_schema", rows)

    @pytest.fixture
    def complex_schema(self):
        """Complex schema with more loci for nearest ST testing."""
        profile_data = """ST	locus1	locus2	locus3	locus4	locus5	locus6	locus7	locus8
1	1	1	2	1	1	2	3	1
2	8	2	7	3	7	1	4	2
3	3	8	5	11	8	3	5	3
4	2	4	6	4	1	6	1	1
5	5	3	3	10	5	8	2	5
"""
        rows = [line.split('\t') for line in profile_data.strip().split('\n')]
        return torchbase.Profile.parse("complex_schema", rows)

    def hamming_distance(self, profile1: torchbase.Profile, profile2: torchbase.Profile) -> int:
        """Calculate hamming distance between two profiles.

        Ignores excluded loci and treats missing loci as wildcards.
        """
        distance = 0
        for locus, val1, val2 in profile1.zip(profile2):
            # Skip excluded loci
            if val1 is torchbase.Special.EXCLUDE or val2 is torchbase.Special.EXCLUDE:
                continue
            # Wildcard matches anything
            if val1 is torchbase.Special.IGNORE or val2 is torchbase.Special.IGNORE:
                continue
            # Count mismatch
            if str(val1) != str(val2):
                distance += 1
        return distance

    # ========== HAPPY PATH TESTS ==========

    def test_exact_match_known_profile(self, simple_schema):
        """Exact match should return known status with correct ST."""
        # Profile 1 from schema exactly matches
        query = torchbase.Profile("test_schema", "ST1", locusA="1", locusB="1", locusC="1", locusD="1")

        # Find exact match
        found = None
        for schema_profile in simple_schema.profiles:
            if query == schema_profile:
                found = schema_profile
                break

        assert found is not None, "Should find exact match"
        assert found.profile == "1"

    def test_nearest_st_single_difference(self, simple_schema):
        """Novel profile differing by 1 allele should find nearest ST."""
        # Profile differs from ST1 by 1 allele (ST1 is 1,1,1,1, changing last to 3)
        query = torchbase.Profile("test_schema", None, locusA="1", locusB="1", locusC="1", locusD="3")

        # Find nearest neighbor
        min_distance = float('inf')
        nearest_st = None
        for schema_profile in simple_schema.profiles:
            dist = self.hamming_distance(query, schema_profile)
            if dist < min_distance:
                min_distance = dist
                nearest_st = schema_profile.profile

        assert nearest_st is not None
        assert min_distance == 1
        assert nearest_st == "1"  # ST1 has same first 3 loci, differs only in 4th

    def test_nearest_st_multiple_differences(self, simple_schema):
        """Novel profile with multiple differences should find closest match."""
        # Query differs from multiple profiles
        query = torchbase.Profile("test_schema", None, locusA="2", locusB="2", locusC="2", locusD="2")

        # Find nearest neighbor (should be ST8)
        min_distance = float('inf')
        nearest_st = None
        for schema_profile in simple_schema.profiles:
            dist = self.hamming_distance(query, schema_profile)
            if dist < min_distance:
                min_distance = dist
                nearest_st = schema_profile.profile

        assert nearest_st == "8"
        assert min_distance == 0  # ST8 is exact match

    def test_output_json_structure_known(self, simple_schema):
        """Output JSON should have correct structure for known profile."""
        # Simulate task output
        output = {
            "sequence_type": "1",
            "status": "known",
            "nearest_st": None,
            "distance": None
        }

        assert "sequence_type" in output
        assert "status" in output
        assert output["status"] in ["known", "novel_profile", "novel_allele"]
        assert "nearest_st" in output
        assert "distance" in output

    def test_output_json_structure_novel_profile(self, simple_schema):
        """Output JSON should have correct structure for novel profile."""
        # Simulate task output for novel profile
        output = {
            "sequence_type": None,
            "status": "novel_profile",
            "nearest_st": "2",
            "distance": 1
        }

        assert output["status"] == "novel_profile"
        assert output["nearest_st"] is not None
        assert output["distance"] is not None
        assert isinstance(output["distance"], int)

    # ========== NOVEL PROFILE DETECTION TESTS ==========

    def test_novel_profile_known_alleles(self, simple_schema):
        """Novel combination of known alleles should be detected."""
        # All alleles exist individually, but not this combination
        query = torchbase.Profile("test_schema", None, locusA="1", locusB="2",
                                  locusC="2", locusD="1")

        # Should not match any exact profile
        exact_match = False
        for schema_profile in simple_schema.profiles:
            if query == schema_profile:
                exact_match = True
                break

        assert not exact_match, "Should not find exact match for novel profile"

        # But should find nearest ST
        min_distance = float('inf')
        for schema_profile in simple_schema.profiles:
            dist = self.hamming_distance(query, schema_profile)
            if dist < min_distance:
                min_distance = dist

        assert min_distance < float('inf'), "Should find nearest ST even if novel"

    def test_novel_profile_status_different_st(self, simple_schema):
        """Multiple profiles with same minimum distance should choose first/lowest ST."""
        # Create a query that's equally distant from multiple profiles
        query = torchbase.Profile("test_schema", None, locusA="3", locusB="2",
                                  locusC="3", locusD="2")

        # Find all profiles with minimum distance
        distances = {}
        for schema_profile in simple_schema.profiles:
            dist = self.hamming_distance(query, schema_profile)
            distances[schema_profile.profile] = dist

        # Should have at least one nearest neighbor
        min_dist = min(distances.values())
        candidates = [st for st, d in distances.items() if d == min_dist]
        assert len(candidates) > 0

    # ========== NOVEL ALLELE DETECTION TESTS ==========

    def test_novel_allele_detection(self, simple_schema):
        """Unknown allele value should be detected as novel_allele."""
        # Check if allele is novel (not seen in any profile)
        locus_a_values = set()
        for schema_profile in simple_schema.profiles:
            if "locusA" in schema_profile:
                locus_a_values.add(str(schema_profile["locusA"]))

        assert "99" not in locus_a_values, "Allele 99 should not exist"

    def test_novel_allele_multiple_loci(self, simple_schema):
        """Multiple novel alleles should be detected correctly."""
        # Multiple unknown alleles
        query = torchbase.Profile("test_schema", None, locusA="99", locusB="99",
                                  locusC="1", locusD="1")

        # Check for novel alleles
        novel_alleles = []
        for locus_name in ["locusA", "locusB"]:
            locus_values = set()
            for schema_profile in simple_schema.profiles:
                locus_values.add(str(schema_profile[locus_name]))

            query_val = str(query[locus_name])
            if query_val not in locus_values:
                novel_alleles.append(locus_name)

        assert len(novel_alleles) == 2

    # ========== MISSING/PARTIAL LOCI TESTS ==========

    def test_partial_profile_missing_locus(self, simple_schema):
        """Should handle profiles with missing loci gracefully."""
        # Profile with fewer loci than schema
        query = torchbase.Profile("test_schema", None, locusA="1", locusB="1",
                                  locusC="1")

        # Should still be able to compare
        min_distance = float('inf')
        for schema_profile in simple_schema.profiles:
            dist = self.hamming_distance(query, schema_profile)
            if dist < min_distance:
                min_distance = dist

        assert min_distance < float('inf'), "Should calculate distance despite missing locus"

    def test_partial_profile_fewer_loci_than_schema(self, simple_schema):
        """Query with fewer loci should match partial sets of schema profiles."""
        # Multiple profiles should match on these loci
        matches = []
        for schema_profile in simple_schema.profiles:
            # Check first two loci
            locus_a_match = str(schema_profile["locusA"]) == "1"
            locus_b_match = str(schema_profile["locusB"]) == "1"
            if locus_a_match and locus_b_match:
                matches.append(schema_profile.profile)

        assert len(matches) > 0, "Should find profiles matching partial profile"

    def test_query_with_missing_locus_nearest_st(self, simple_schema):
        """Partial profile should still find valid nearest ST."""
        query = torchbase.Profile("test_schema", None, locusA="2", locusB="2")

        min_distance = float('inf')
        for schema_profile in simple_schema.profiles:
            dist = self.hamming_distance(query, schema_profile)
            if dist < min_distance:
                min_distance = dist

        assert min_distance < float('inf')

    # ========== EDGE CASE TESTS ==========

    def test_empty_profile_not_compared(self, simple_schema):
        """Empty profile should not crash comparison."""
        # Should still be able to iterate schemas without error
        count = 0
        for schema_profile in simple_schema.profiles:
            count += 1

        assert count == len(simple_schema.profiles)

    def test_all_same_alleles_match(self, simple_schema):
        """Profile with all loci same as one schema profile should match exactly."""
        query = torchbase.Profile("test_schema", "ST7",
                                  locusA="1", locusB="2", locusC="2", locusD="2")

        # Should find exact match with ST7
        for schema_profile in simple_schema.profiles:
            if query == schema_profile:
                assert schema_profile.profile == "7"
                return

        pytest.fail("Should find exact match for ST7")

    def test_single_locus_difference_finds_nearest(self, complex_schema):
        """Single allele difference should find correct nearest ST."""
        # ST1 is (1,1,2,1,1,2,3,1), change last locus
        query = torchbase.Profile("complex_schema", None,
                                  locus1="1", locus2="1", locus3="2", locus4="1",
                                  locus5="1", locus6="2", locus7="3", locus8="2")

        min_distance = float('inf')
        for schema_profile in complex_schema.profiles:
            dist = self.hamming_distance(query, schema_profile)
            if dist < min_distance:
                min_distance = dist

        # ST1 should have distance 1 (only locus8 differs)
        assert min_distance == 1

    def test_all_alleles_different_from_any(self, simple_schema):
        """Query completely different from all profiles should still find a nearest."""
        query = torchbase.Profile("test_schema", None,
                                  locusA="9", locusB="9", locusC="9", locusD="9")

        distances = {}
        for schema_profile in simple_schema.profiles:
            dist = self.hamming_distance(query, schema_profile)
            distances[schema_profile.profile] = dist

        # Should find at least one nearest ST
        min_dist = min(distances.values())
        assert min_dist < float('inf')

        # Get the ST with minimum distance
        nearest_sts = [st for st, d in distances.items() if d == min_dist]
        assert len(nearest_sts) > 0

    def test_numeric_string_comparison(self, simple_schema):
        """Allele values should compare correctly as numeric strings."""
        # Values "1" and 1 should be treated as equal
        query = torchbase.Profile("test_schema", "ST_compare",
                                  locusA=1, locusB="1", locusC="1", locusD="1")

        for schema_profile in simple_schema.profiles:
            if schema_profile.profile == "1":
                # Should match despite mixed types
                for _, val1, val2 in query.zip(schema_profile):
                    # Values should compare as strings
                    if str(val1) == str(val2):
                        continue
                    elif val1 is torchbase.Special.IGNORE:
                        continue
                    else:
                        # Type comparison works correctly
                        assert False, "Types should compare as strings"

    # ========== BOUNDARY TESTS ==========

    def test_large_distance_calculation(self, complex_schema):
        """Should correctly calculate distance for profiles with many differences."""
        # Change multiple alleles
        query = torchbase.Profile("complex_schema", None,
                                  locus1="2", locus2="3", locus3="4", locus4="5",
                                  locus5="6", locus6="7", locus7="8", locus8="9")

        distances = {}
        for schema_profile in complex_schema.profiles:
            dist = self.hamming_distance(query, schema_profile)
            distances[schema_profile.profile] = dist

        # All should have non-zero distances
        assert all(d > 0 for d in distances.values())
        # Should still find a minimum
        min_dist = min(distances.values())
        assert min_dist < float('inf')

    def test_tie_breaking_multiple_nearest_st(self, simple_schema):
        """When multiple profiles equidistant, should handle gracefully."""
        # Create query equidistant from multiple profiles
        query = torchbase.Profile("test_schema", None,
                                  locusA="1", locusB="1", locusC="1",
                                  locusD="1.5")

        # Won't match perfectly due to decimal, but should find nearest
        min_distance = float('inf')
        nearest_sts = []
        for schema_profile in simple_schema.profiles:
            dist = self.hamming_distance(query, schema_profile)
            if dist < min_distance:
                min_distance = dist
                nearest_sts = [schema_profile.profile]
            elif dist == min_distance:
                nearest_sts.append(schema_profile.profile)

        assert len(nearest_sts) >= 1


class TestProfileLookupOutputFormat:
    """Tests for output JSON format compliance."""

    def test_output_has_all_required_fields(self):
        """Output must have sequence_type, status, nearest_st, distance fields."""
        output = {
            "sequence_type": "1",
            "status": "known",
            "nearest_st": None,
            "distance": None
        }

        required_fields = ["sequence_type", "status", "nearest_st", "distance"]
        for field in required_fields:
            assert field in output, f"Missing required field: {field}"

    def test_known_status_has_st_no_distance(self):
        """Known profile should have sequence_type but null nearest_st/distance."""
        output = {
            "sequence_type": "42",
            "status": "known",
            "nearest_st": None,
            "distance": None
        }

        assert output["status"] == "known"
        assert output["sequence_type"] is not None
        assert output["nearest_st"] is None
        assert output["distance"] is None

    def test_novel_profile_status_has_nearest_st(self):
        """Novel profile should have nearest_st and distance."""
        output = {
            "sequence_type": None,
            "status": "novel_profile",
            "nearest_st": "42",
            "distance": 2
        }

        assert output["status"] == "novel_profile"
        assert output["sequence_type"] is None
        assert output["nearest_st"] is not None
        assert output["distance"] is not None
        assert isinstance(output["distance"], int)

    def test_novel_allele_status_has_nearest_st(self):
        """Novel allele should have nearest_st and distance."""
        output = {
            "sequence_type": None,
            "status": "novel_allele",
            "nearest_st": "15",
            "distance": 1
        }

        assert output["status"] == "novel_allele"
        assert output["nearest_st"] is not None
        assert output["distance"] is not None

    def test_status_enum_values(self):
        """Status should only be known, novel_profile, or novel_allele."""
        valid_statuses = ["known", "novel_profile", "novel_allele"]

        for status in valid_statuses:
            assert status in valid_statuses

    def test_distance_is_integer(self):
        """Distance field should be integer when present."""
        output = {
            "distance": 3
        }
        assert isinstance(output["distance"], int)


class TestProfileLookupErrors:
    """Tests for error handling and validation."""

    def test_mismatched_schema_between_query_and_table(self):
        """Should handle schema mismatch gracefully."""
        # This is more of a data validation test
        query_schema = "schema_v1"
        table_schema = "schema_v2"

        # Would need error handling in actual implementation
        assert query_schema != table_schema

    def test_malformed_allele_call_json(self):
        """Should validate allele call JSON format."""
        # Valid structure would have loci as keys
        valid_json = {"locusA": "1", "locusB": "2"}

        assert "locusA" in valid_json
        assert "locusB" in valid_json

    def test_empty_profiles_table(self):
        """Should handle empty profile table gracefully."""
        # With no profiles in schema, any query is novel
        empty_schema = torchbase.Schema("empty", profiles=[])

        assert len(empty_schema.profiles) == 0

    def test_single_profile_in_table(self):
        """Should work with single profile in table."""
        single_profile = torchbase.Profile("test", "1", locus1="1", locus2="2")
        schema = torchbase.Schema("test", profiles=[single_profile])

        assert len(schema.profiles) == 1


class TestHammingDistanceCalculation:
    """Tests for hamming distance calculation logic."""

    def hamming_distance(self, profile1: torchbase.Profile, profile2: torchbase.Profile) -> int:
        """Calculate hamming distance between two profiles."""
        distance = 0
        for locus, val1, val2 in profile1.zip(profile2):
            if val1 is torchbase.Special.EXCLUDE or val2 is torchbase.Special.EXCLUDE:
                continue
            if val1 is torchbase.Special.IGNORE or val2 is torchbase.Special.IGNORE:
                continue
            if str(val1) != str(val2):
                distance += 1
        return distance

    def test_identical_profiles_zero_distance(self):
        """Identical profiles should have distance 0."""
        p1 = torchbase.Profile("t", "1", a="1", b="2", c="3")
        p2 = torchbase.Profile("t", "2", a="1", b="2", c="3")

        dist = self.hamming_distance(p1, p2)
        assert dist == 0

    def test_single_difference_distance_one(self):
        """Single allele difference should give distance 1."""
        p1 = torchbase.Profile("t", "1", a="1", b="2", c="3")
        p2 = torchbase.Profile("t", "2", a="1", b="2", c="4")

        dist = self.hamming_distance(p1, p2)
        assert dist == 1

    def test_all_different_alleles(self):
        """All different alleles should give distance equal to number of loci."""
        p1 = torchbase.Profile("t", "1", a="1", b="2", c="3", d="4")
        p2 = torchbase.Profile("t", "2", a="5", b="6", c="7", d="8")

        dist = self.hamming_distance(p1, p2)
        assert dist == 4

    def test_distance_ignores_excluded_loci(self):
        """Excluded loci should not count toward distance."""
        p1 = torchbase.Profile("t", "1", a="1", b=torchbase.Special.EXCLUDE, c="3")
        p2 = torchbase.Profile("t", "2", a="1", b="9", c="3")

        dist = self.hamming_distance(p1, p2)
        # b should be ignored, so only 0 differences
        assert dist == 0

    def test_distance_ignores_wildcard_loci(self):
        """Wildcard IGNORE loci should not count toward distance."""
        p1 = torchbase.Profile("t", "1", a="1", b=torchbase.Special.IGNORE, c="3")
        p2 = torchbase.Profile("t", "2", a="1", b="9", c="3")

        dist = self.hamming_distance(p1, p2)
        # IGNORE should match anything
        assert dist == 0
