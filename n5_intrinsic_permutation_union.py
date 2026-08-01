#!/usr/bin/env python3
"""Permutation helpers shared by the five-player family searches."""

from __future__ import annotations


def permute_mask(mask, permutation):
    return sum(
        1 << permutation[player]
        for player in range(len(permutation))
        if mask >> player & 1
    )


def permute_game(game, permutation):
    image = [None] * len(game)
    for coalition, value in enumerate(game):
        image[permute_mask(coalition, permutation)] = value
    return tuple(image)
