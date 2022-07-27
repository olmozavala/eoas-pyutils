from proc_utils.comp_fields import coriolis
import pytest

def test_coriolis():
    '''
    Tests the coriolis function
    '''
    assert coriolis(0) == 0
    assert coriolis(45) == pytest.approx(1.0284e-04, rel=1e-3)
    assert coriolis([0, -10]) == pytest.approx([0, -2.5256e-5], rel=1e-3)