# Copyright (C) 2024 Federico Battista, Fabrizio Rossi, Stefano Smriglio
# SPDX-License-Identifier: GPL-3.0-or-later
# Tests for graph generation, I/O utilities, and clique cover algorithms.

import os
import tempfile
import itertools
import networkx as nx
import pytest

from pyModules.Graphs import (
    web, antiweb, wheel, cycle, random_graph,
    greedy_clique_cover, greedy_clique_cover_letchford_et_al,
    write_graph_to_dimacs, read_graph_from_dimacs,
)

# ---------------------------------------------------------------------------
# DIMACS I/O round-trip
# ---------------------------------------------------------------------------

class TestDIMACSIO:

    def _roundtrip(self, G):
        with tempfile.NamedTemporaryFile(suffix='.stb', delete=False) as f:
            path = f.name
        try:
            write_graph_to_dimacs(G, path)
            G2 = read_graph_from_dimacs(path)
        finally:
            os.unlink(path)
        return G2

    def test_roundtrip_cycle(self):
        G = cycle(6)
        G2 = self._roundtrip(G)
        assert G.number_of_nodes() == G2.number_of_nodes()
        assert G.number_of_edges() == G2.number_of_edges()

    def test_roundtrip_web(self):
        G = web(7, 2)
        G2 = self._roundtrip(G)
        assert G.number_of_edges() == G2.number_of_edges()

    def test_dimacs_format(self):
        G = cycle(4)
        with tempfile.NamedTemporaryFile(suffix='.stb', delete=False, mode='w') as f:
            path = f.name
        try:
            write_graph_to_dimacs(G, path)
            with open(path) as f:
                header = f.readline()
            assert header.startswith('p edge')
        finally:
            os.unlink(path)


# ---------------------------------------------------------------------------
# Clique cover algorithms
# ---------------------------------------------------------------------------

class TestCliqueCover:

    def _all_edges_covered(self, G, cover):
        """Check that every edge of G is covered by at least one clique."""
        for u, v in G.edges():
            covered = any(u in clique and v in clique for clique in cover)
            assert covered, "Edge (%d, %d) not covered" % (u, v)

    def _each_clique_is_clique(self, G, cover):
        """Check that each element of the cover is a clique in G."""
        for clique in cover:
            for u, v in itertools.combinations(clique, 2):
                assert G.has_edge(u, v), "(%d, %d) not an edge but in clique" % (u, v)

    def test_greedy_cover_covers_all_edges(self):
        G = web(7, 2)
        cover = greedy_clique_cover(G)
        self._all_edges_covered(G, cover)

    def test_greedy_cover_each_set_is_clique(self):
        G = web(7, 2)
        cover = greedy_clique_cover(G)
        self._each_clique_is_clique(G, cover)

    def test_letchford_cover_covers_all_edges(self):
        G = web(7, 2)
        cover = greedy_clique_cover_letchford_et_al(G)
        self._all_edges_covered(G, cover)

    def test_letchford_cover_each_set_is_clique(self):
        G = web(7, 2)
        cover = greedy_clique_cover_letchford_et_al(G)
        self._each_clique_is_clique(G, cover)

    def test_empty_graph_returns_empty_cover(self):
        G = nx.Graph()
        G.add_nodes_from(range(4))
        cover = greedy_clique_cover(G)
        assert cover == []

    def test_complete_graph_covered_by_one_clique(self):
        G = nx.complete_graph(5)
        cover = greedy_clique_cover_letchford_et_al(G)
        assert len(cover) == 1
        assert set(cover[0]) == set(range(5))
