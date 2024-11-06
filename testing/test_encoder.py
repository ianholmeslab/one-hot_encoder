import pytest
from src.encoder import Encoder

class TestEncoder():

    x = 17

    def test_tautology(self):
        assert self.x == 17
        assert str(self.x) != "test_string"

    def test_encoder_accession(self):
        assert self.x % 5 != 0
        assert Encoder.pytest_test_method()
