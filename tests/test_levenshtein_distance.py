#!/usr/bin/env python
# -*- coding: utf-8 -*-

import Levenshtein
import unittest

def test_empty_string():
    """
    when both strings are empty this is a perfect match
    """
    assert Levenshtein.distance(b"", b"") == 0
    assert Levenshtein.distance(u"", u"") == 0

def test_simple():
    """
    some very simple tests using supported string types
    bytes/unicode
    to catch relatively obvious implementation errors
    """
    assert Levenshtein.distance(b"ABCD", b"AF") == 3
    assert Levenshtein.distance(u"ABCD", u"AF") == 3
    assert Levenshtein.distance(b"ABCD", b"ABCD") == 0
    assert Levenshtein.distance(u"ABCD", u"ABCD") == 0

def test_simple_unicode_tests():
    """
    some very simple tests using unicode
    to catch relatively obvious implementation errors
    """
    s1 = u"ÁÄ"
    s2 = u"ABCD"
    assert Levenshtein.distance(s1, s2) == 4
    assert Levenshtein.distance(s1, s1) == 0
