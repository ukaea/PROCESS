import pytest

profile = NeutronFluxProfile(
    1.0, #flux
    "dummy-FW-mat",
    0.01, #1cm
    "dummy-BZ-mat",
    0.11, #10cm thick blanket
    n_groups=2,
    )
def test_has_fluxes():
    assert hasattr(profile, "neutron_flux_at")
    assert hasattr(profile, "groupwise_neutron_flux_at")
    assert hasattr(profile, "neutron_flux_fw")
    assert hasattr(profile, "groupwise_neutron_flux_fw")
    assert hasattr(profile, "neutron_flux_bz")
    assert hasattr(profile, "groupwise_neutron_flux_bz")

def test_reactions():
    assert hasattr(profile, "groupwise_reaction_rate_fw")
    assert hasattr(profile, "reaction_rate_fw")
    assert hasattr(profile, "groupwise_reaction_rate_bz")
    assert hasattr(profile, "reaction_rate_bz")
    assert hasattr(profile, "flux_fw2bz")
    assert hasattr(profile, "flux_escaped")