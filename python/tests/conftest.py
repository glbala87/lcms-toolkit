"""Pytest configuration and shared fixtures."""

import pytest
import numpy as np


@pytest.fixture
def random_seed():
    """Set a fixed random seed for reproducible tests."""
    np.random.seed(42)
    yield
    np.random.seed(None)
