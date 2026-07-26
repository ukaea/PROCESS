import pytest

from process.models.physics.l_h_transition import PlasmaConfinementTransition


@pytest.fixture
def physics(process_models):
    """Provides Physics object for testing.

    :returns: initialised Physics object
    :rtype: process.physics.Physics
    """
    return process_models.physics


@pytest.mark.parametrize(
    ("a", "b", "c", "expected"),
    [
        (1.0, 5.0, 6.2, 86.49000000000001),
    ],
)
def test_calculate_iter1996_nominal(a, b, c, expected):
    assert PlasmaConfinementTransition().calculate_iter1996_nominal(
        a, b, c
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "c", "expected"),
    [
        (1.0, 5.0, 6.2, 189.5394231299933),
    ],
)
def test_calculate_iter1996_upper(a, b, c, expected):
    assert PlasmaConfinementTransition().calculate_iter1996_upper(
        a, b, c
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "c", "expected"),
    [
        (1.0, 5.0, 6.2, 39.46682952353113),
    ],
)
def test_calculate_iter1996_lower(a, b, c, expected):
    assert PlasmaConfinementTransition().calculate_iter1996_lower(
        a, b, c
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "c", "expected"),
    [
        (1.0, 5.0, 6.2, 131.12022710987677),
    ],
)
def test_calculate_snipes1997_iter(a, b, c, expected):
    assert PlasmaConfinementTransition().calculate_snipes1997_iter(
        a, b, c
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "c", "d", "expected"),
    [
        (1.0, 5.0, 6.2, 1.8, 105.48455140282398),
    ],
)
def test_calculate_snipes1997_kappa(a, b, c, d, expected):
    assert PlasmaConfinementTransition().calculate_snipes1997_kappa(
        a, b, c, d
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "c", "d", "expected"),
    [
        (1.0, 5.0, 100.0, 2.0, 13.542309512892546),
    ],
)
def test_calculate_martin08_nominal(a, b, c, d, expected):
    assert PlasmaConfinementTransition().calculate_martin08_nominal(
        a, b, c, d
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "c", "d", "expected"),
    [
        (1.0, 5.0, 100.0, 2.0, 16.474587948461522),
    ],
)
def test_calculate_martin08_upper(a, b, c, d, expected):
    assert PlasmaConfinementTransition().calculate_martin08_upper(
        a, b, c, d
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "c", "d", "expected"),
    [
        (1.0, 5.0, 100.0, 2.0, 11.131941359980322),
    ],
)
def test_calculate_martin08_lower(a, b, c, d, expected):
    assert PlasmaConfinementTransition().calculate_martin08_lower(
        a, b, c, d
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "c", "d", "e", "expected"),
    [
        (1.0, 5.0, 6.2, 2.0, 2.0, 57.765659228267175),
    ],
)
def test_calculate_snipes2000_nominal(a, b, c, d, e, expected):
    assert PlasmaConfinementTransition().calculate_snipes2000_nominal(
        a, b, c, d, e
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "c", "d", "e", "expected"),
    [
        (1.0, 5.0, 6.2, 2.0, 2.0, 81.45741303978973),
    ],
)
def test_calculate_snipes2000_upper(a, b, c, d, e, expected):
    assert PlasmaConfinementTransition().calculate_snipes2000_upper(
        a, b, c, d, e
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "c", "d", "e", "expected"),
    [
        (1.0, 5.0, 6.2, 2.0, 2.0, 40.636940606933976),
    ],
)
def test_calculate_snipes2000_lower(a, b, c, d, e, expected):
    assert PlasmaConfinementTransition().calculate_snipes2000_lower(
        a, b, c, d, e
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "c", "d", "expected"),
    [
        (1.0, 5.0, 6.2, 2.0, 29.515866708312462),
    ],
)
def test_calculate_snipes2000_closed_divertor_nominal(a, b, c, d, expected):
    assert PlasmaConfinementTransition().calculate_snipes2000_closed_divertor_nominal(
        a, b, c, d
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "c", "d", "expected"),
    [
        (1.0, 5.0, 6.2, 2.0, 40.414673189823354),
    ],
)
def test_calculate_snipes2000_closed_divertor_upper(a, b, c, d, expected):
    assert PlasmaConfinementTransition().calculate_snipes2000_closed_divertor_upper(
        a, b, c, d
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "c", "d", "expected"),
    [
        (1.0, 5.0, 6.2, 2.0, 21.404993867161203),
    ],
)
def test_calculate_snipes2000_closed_divertor_lower(a, b, c, d, expected):
    assert PlasmaConfinementTransition().calculate_snipes2000_closed_divertor_lower(
        a, b, c, d
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "expected"),
    [
        (1e6, 1.0, 2.11),
    ],
)
def test_calculate_hubbard2012_nominal(a, b, expected):
    assert PlasmaConfinementTransition().calculate_hubbard2012_nominal(
        a, b
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "expected"),
    [
        (1e6, 1.0, 2.11),
    ],
)
def test_calculate_hubbard2012_upper(a, b, expected):
    assert PlasmaConfinementTransition().calculate_hubbard2012_upper(
        a, b
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "expected"),
    [
        (1e6, 1.0, 2.11),
    ],
)
def test_calculate_hubbard2012_lower(a, b, expected):
    assert PlasmaConfinementTransition().calculate_hubbard2012_lower(
        a, b
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "c", "expected"),
    [
        (1.0, 100.0, 5.0, 24.61768530477778),
    ],
)
def test_calculate_hubbard2017(a, b, c, expected):
    assert PlasmaConfinementTransition().calculate_hubbard2017(a, b, c) == pytest.approx(
        expected
    )


@pytest.mark.parametrize(
    ("a", "b", "c", "d", "e", "expected"),
    [
        (1.0, 5.0, 100.0, 2.0, 2.5, 13.593852185804629),
    ],
)
def test_calculate_martin08_aspect_nominal(a, b, c, d, e, expected):
    assert PlasmaConfinementTransition().calculate_martin08_aspect_nominal(
        a, b, c, d, e
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "c", "d", "e", "expected"),
    [
        (1.0, 5.0, 100.0, 2.0, 2.5, 16.537291012306024),
    ],
)
def test_calculate_martin08_aspect_upper(a, b, c, d, e, expected):
    assert PlasmaConfinementTransition().calculate_martin08_aspect_upper(
        a, b, c, d, e
    ) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("a", "b", "c", "d", "e", "expected"),
    [
        (1.0, 5.0, 100.0, 2.0, 2.5, 11.174310057272887),
    ],
)
def test_calculate_martin08_aspect_lower(a, b, c, d, e, expected):
    assert PlasmaConfinementTransition().calculate_martin08_aspect_lower(
        a, b, c, d, e
    ) == pytest.approx(expected)
