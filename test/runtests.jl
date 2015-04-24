using GibbsSeaWater
using Compat
using Base.Test

testTolerances = @compat Dict(:abs_pressure_from_p => 2.300031483173370e-004,
                              :adiabatic_lapse_rate_from_ct => 5.699743845487985e-019,
                              :adiabatic_lapse_rate_from_t => 5.699738675609156e-019,
                              :adiabatic_lapse_rate_t_exact => 5.699738675609156e-019,
                              :alpha => 8.264713918662917e-015,
                              :alpha_ct_exact => 8.257102480598889e-015,
                              :alpha_ct_exact_rab => 8.257102480598889e-015,
                              :alpha_on_beta => 1.052907760978883e-11,
                              :alpha_rab => 8.264713918662917e-015,
                              :alpha_wrt_ct_t_exact => 8.257094010269417e-015,
                              :alpha_wrt_pt_t_exact => 8.599068316130290e-015,
                              :alpha_wrt_t_exact => 8.594856868316542e-015,
                              :beta => 1.846179459308317e-015,
                              :beta_const_ct_t_exact => 1.847372081698051e-015,
                              :beta_const_pt_t_exact => 1.805738718274608e-015,
                              :beta_const_t_exact => 1.804871356536619e-015,
                              :beta_ct_exact => 1.847372081698051e-015,
                              :beta_ct_exact_rab => 1.847372081698051e-015,
                              :beta_rab => 1.846179459308317e-015,
                              :brinesa_ct => 1.300080043620255e-010,
                              :brinesa_t => 1.300293206440983e-010,
                              :cabbeling => 1.634722766223964e-16,
                              :c_from_sp => 6.163816124171717e-010,
                              :chem_potential_salt_t_exact => 4.805315256817266e-007,
                              :chem_potential_t_exact => 1.317864928296331e-009,
                              :chem_potential_water_t_exact => 4.811790859093890e-007,
                              :cp_t_exact => 2.813976607285440e-009,
                              :ct_freezing => 2.257127817983928e-011,
                              :ct_from_entropy => 6.261107188265669e-010,
                              :ct_from_pt => 6.261107188265669e-010,
                              :ct_from_rho => 6.298552790440226e-010,
                              :ct_from_rho_exact => 6.280096442878858e-010,
                              :ct_from_t => 6.261107188265669e-010,
                              :ct_maxdensity => 6.688338771709823e-011,
                              :ct_maxdensity_exact => 5.839062566792563e-011,
                              :ct_pt => 2.964295475749168e-013,
                              :ct_pt_pt => 5.551115123125783e-014,
                              :ct_sa => 1.006231122729906e-012,
                              :ct_sa_pt => 1.457167719820518e-014,
                              :ct_sa_sa => 1.431146867680866e-014,
                              :deltasa_from_sp => 6.963318810448982e-013,
                              :deltasa_atlas => 6.945514042372425e-013,
                              :depth_from_z => 2.287223921371151e-008,
                              :distance => 4.470348358154297e-008,
                              :drho_dct => 8.315403920988729e-12,
                              :drho_dp => 1.782157321023048e-18,
                              :drho_dsa => 1.943112337698949e-12,
                              :dynamic_enthalpy => 2.288754734930485e-007,
                              :dynamic_enthalpy_ct_exact => 2.288797986693680e-007,
                              :dynamic_enthalpy_t_exact => 2.288943505845964e-007,
                              :enthalpy => 2.499356924090534e-006,
                              :enthalpy_ct_exact => 2.499349648132920e-006,
                              :enthalpy_diff => 1.154347728515859e-010,
                              :enthalpy_diff_ct_exact => 7.275957614183426e-011,
                              :enthalpy_t_exact => 2.499349648132920e-006,
                              :entropy_from_ct => 9.028163105995191e-009,
                              :entropy_from_pt => 9.028163105995191e-009,
                              :entropy_from_t => 9.028163105995191e-009,
                              :entropy_t_exact => 9.028163105995191e-009,
                              :eta_ct => 3.137579085432662e-011,
                              :eta_ct_ct => 2.381358998881922e-013,
                              :eta_sa => 4.653527563291959e-012,
                              :eta_sa_ct => 4.981922535098049e-014,
                              :eta_sa_sa => 6.995931611797346e-013,
                              :f => 1.000000000000000e-015,
                              :fdelta => 2.702939055302528e-014,
                              :g_ct => 2.854610769986721e-009,
                              :geo_strf_cunningham => 9.926097845891491e-007,
                              :geo_strf_dyn_height => 9.969444363377988e-007,
                              :geo_strf_dyn_height_pc => 4.210555459849275e-008,
                              :geo_strf_dyn_height_pc_p_mid => 1.000000000000000e-015,
                              :geo_strf_isopycnal => 4.297791695861974e-007,
                              :geo_strf_isopycnal_pc => 8.540336438045415e-009,
                              :geo_strf_isopycnal_pc_p_mid => 1.000000000000000e-015,
                              :geo_strf_montgomery => 9.915690188933013e-007,
                              :geo_strf_velocity => 8.024344404222727e-009,
                              :geo_strf_velocity_mid_lat => 1.000000000000000e-015,
                              :geo_strf_velocity_mid_long => 1.000000000000000e-015,
                              :grav => 5.329070518200751e-014,
                              :h_ct => 3.160494088660926e-010,
                              :h_ct_ct => 1.135647131889073e-011,
                              :helmholtz_energy_t_exact => 2.440137905068696e-007,
                              :h_p => 2.818925648462312e-016,
                              :h_sa => 2.371365326325758e-010,
                              :h_sa_ct => 2.188554892867956e-012,
                              :h_sa_sa => 7.264189250122399e-013,
                              :internal_energy => 2.499342372175306e-006,
                              :internal_energy_ct_exact => 2.499335096217692e-006,
                              :internal_energy_t_exact => 2.499335096217692e-006,
                              :ionic_strength_from_sa => 2.768674178810215e-012,
                              :ipvfn2 => 3.474816878679121e-009,
                              :isochoric_heat_cap_t_exact => 1.614353095646948e-009,
                              :isopycnal_slope_ratio => 3.781384094736495e-010,
                              :kappa => 1.717743939542931e-21,
                              :kappa_const_t_exact => 1.697064424229105e-021,
                              :kappa_t_exact => 1.712677458291044e-021,
                              :latentheat_evap_ct => 1.455657184123993e-006,
                              :latentheat_evap_t => 1.443084329366684e-006,
                              :latentheat_melting => 6.286427378654480e-008,
                              :molality_from_sa => 4.446665258228677e-012,
                              :n2 => 1.578186366313350e-014,
                              :ntpptct => 5.024869409453459e-013,
                              :osmotic_coefficient_t_exact => 3.583799923490005e-013,
                              :osmotic_pressure_t_exact => 1.023465756588848e-009,
                              :p_from_abs_pressure => 2.300066626048647e-008,
                              :p_from_z => 2.301931090187281e-008,
                              :p_mid_g_ct => 2.300021151313558e-008,
                              :p_mid_ipvfn2 => 2.300021151313558e-008,
                              :p_mid_n2 => 2.300021151313558e-008,
                              :p_mid_tursr => 2.300021151313558e-008,
                              :pot_enthalpy_from_pt => 2.499356924090534e-006,
                              :pot_rho_t_exact => 2.929709808086045e-010,
                              :pt0_from_t => 6.054037271496782e-010,
                              :pt_ct => 2.733369086627135e-013,
                              :pt_ct_ct => 4.440892098500626e-014,
                              :pt_from_ct => 6.054037271496782e-010,
                              :pt_from_entropy => 6.054072798633570e-010,
                              :pt_from_t => 6.054037271496782e-010,
                              :pt_sa => 9.670458878119348e-013,
                              :pt_sa_ct => 1.199127602768968e-014,
                              :pt_sa_sa => 2.081668171172169e-014,
                              :r_from_sp => 1.436317731418058e-011,
                              :rho => 2.945625965367071e-010,
                              :rho_ct_exact => 2.944489096989855e-010,
                              :rho_ct_exact_rab => 2.944489096989855e-010,
                              :rho_rab => 2.944489096989855e-010,
                              :rho_t_exact => 2.944489096989855e-010,
                              :rsubrho => 1.709803143512545e-008,
                              :sa_from_rho => 1.308677610722953e-010,
                              :sa_from_rho_ct_exact => 1.306119656874216e-010,
                              :sa_from_rho_t_exact => 1.304769625676272e-010,
                              :sa_from_sp => 1.300080043620255e-010,
                              :sa_from_sp_baltic => 1.300080043620255e-010,
                              :sa_from_sstar => 1.300222152167407e-010,
                              :sa_sa_sstar_from_sp => 1.300008989346679e-010,
                              :sigma0 => 2.933120413217694e-010,
                              :sigma0_ct_exact => 2.929709808086045e-010,
                              :sigma0_pt0_exact => 2.929709808086045e-010,
                              :sigma1 => 2.999058779096231e-010,
                              :sigma1_ct_exact => 2.994511305587366e-010,
                              :sigma2 => 3.060449671465904e-010,
                              :sigma2_ct_exact => 3.055902197957039e-010,
                              :sigma3 => 3.119566827081144e-010,
                              :sigma3_ct_exact => 3.117293090326712e-010,
                              :sigma4 => 3.180957719450817e-010,
                              :sigma4_ct_exact => 3.176410245941952e-010,
                              :sound_speed => 2.596152626210824e-009,
                              :sound_speed_ct_exact => 2.590240910649300e-009,
                              :sound_speed_t_exact => 2.590240910649300e-009,
                              :specvol_anom => 2.810252031082428e-016,
                              :specvol_anom_ct_exact => 2.805915222392486e-016,
                              :specvol_anom_t_exact => 2.805915222392486e-016,
                              :specvol => 2.821094052807283e-016,
                              :specvol_ct_exact => 2.818925648462312e-016,
                              :specvol_t_exact => 2.818925648462312e-016,
                              :sp_from_c => 1.297193463756230e-010,
                              :sp_from_r => 2.681299626772216e-012,
                              :sp_from_sa => 1.297113527698457e-010,
                              :sp_from_sa_baltic => 1.297113527698457e-010,
                              :sp_from_sr => 1.297122409482654e-010,
                              :sp_from_sstar => 1.297122409482654e-010,
                              :sp_salinometer => 1.297131291266851e-010,
                              :sr_from_sp => 1.303233077010191e-010,
                              :sstar_from_sa => 1.300008989346679e-010,
                              :sstar_from_sp => 1.300008989346679e-010,
                              :sstar_sa_sstar_from_sp => 1.300008989346679e-010,
                              :steric_height => 1.017674460257467e-007,
                              :t90_from_t48 => 5.997407015456702e-010,
                              :t90_from_t68 => 5.998579410970706e-010,
                              :t_freezing => 2.157829470661454e-011,
                              :t_from_ct => 6.000142604989378e-010,
                              :t_from_rho_exact => 6.032974120273593e-010,
                              :thermobaric => 4.890907320111513e-23,
                              :t_maxdensity_exact => 6.274447628129565e-011,
                              :tu => 2.190718007000214e-008,
                              :z_from_depth => 2.287223921371151e-008,
                              :z_from_p => 2.287223921371151e-008)

# Test data
sp =  35.5e0
sa = 35.7e0
sstar = 35.5e0
sr = 35.5e0
t = 15e0
ct = 20e0
pt = 15e0
p = 300e0
p_bs = 50e0
p_ref = 100e0
c = 43.6e0
lon = 260e0
long_bs = 20e0
lat = 16e0
lat_bs = 60e0
saturation_fraction = 0.5e0

function run_gsw_test(funcsym, expected, args...)
    tol = testTolerances[funcsym]
    res = eval(funcsym)(args...)
    @test abs(res-expected) < tol
end

# Tests copied from gsw_check_functions.c

# Practical Salinity, PSS-78
run_gsw_test(:sp_from_c, 35.500961780774482e0, c,t,p);
run_gsw_test(:c_from_sp, 43.598945605280484e0, sp,t,p);

# Absolute Salinity, Preformed Salinity and Conservative Temperature
run_gsw_test(:sa_from_sp, 35.671358392019094e0, sp,p,lon,lat);
run_gsw_test(:sstar_from_sp, 35.866946753006239e0, sa,p,lon,lat);
run_gsw_test(:ct_from_t,  14.930280459895560e0, sa, t, p);

# Other conversions between temperatures, salinities, entropy, pressure and height
run_gsw_test(:deltasa_from_sp, 3.96067773336028495e-3, sp,p,lon,lat);
run_gsw_test(:sr_from_sp, 35.667397714285734e0, sp);
run_gsw_test(:sp_from_sr, 35.333387933015295e0, sr);
run_gsw_test(:sp_from_sa, 35.528504019167094e0, sa,p,lon,lat);
run_gsw_test(:sstar_from_sa, 35.694648791860907e0, sa,p,lon,lat);
run_gsw_test(:sp_from_sstar, 35.334761242083573e0, sstar,p,lon,lat);
run_gsw_test(:sa_from_sstar, 35.505322027120805e0, sstar,p,lon,lat);
run_gsw_test(:pt_from_ct, 20.023899375975017e0, sa,ct);
run_gsw_test(:t_from_ct, 20.079820359223014e0, sa,ct,p);
run_gsw_test(:ct_from_pt, 14.976021403957613e0, sa,pt);
run_gsw_test(:pt0_from_t, 14.954241363902305e0, sa,t,p);
run_gsw_test(:pt_from_t, 14.969381237883740e0, sa,t,p,p_ref);
run_gsw_test(:z_from_p, -2.980161553316402e2, p,lat);
run_gsw_test(:entropy_from_t, 212.30166821093002e0, sa,t,p);
run_gsw_test(:adiabatic_lapse_rate_from_ct, 1.877941744191155e-8, sa,ct,p);

# Density and enthalpy, based on the 48-term expression for density
run_gsw_test(:rho, 1026.4562376198473e0, sa,ct,p);
run_gsw_test(:alpha, 2.62460550806784356e-4, sa,ct,p);
run_gsw_test(:beta, 7.29314455934463365e-4, sa,ct,p);
run_gsw_test(:alpha_on_beta, 0.359872958325632e0, sa,ct,p);

let
    # WIP: rho_first_derivatives expects a pointer, so need to pass an array for now
    drho_dsa = Float64[0.0]
    drho_dct = Float64[0.0]
    drho_dp = Float64[0.0]

    rho_first_derivatives(sa, ct, p, drho_dsa, drho_dct, drho_dp);
    drho_dsa_error = abs(drho_dsa[1] - 0.748609372480258e0);
    drho_dct_error = abs(drho_dct[1] + 0.269404269504765e0);
    drho_dp_error = abs(drho_dp[1] - 4.287533235942749e-7);
    @test drho_dsa_error < testTolerances[:drho_dsa]
    @test drho_dct_error < testTolerances[:drho_dct]
    @test drho_dp_error < testTolerances[:drho_dp]
end

run_gsw_test(:specvol, 9.74225654586897711e-4, sa,ct,p);
run_gsw_test(:specvol_anom, 2.90948181201264571e-6, sa,ct,p);
run_gsw_test(:sigma0, 25.165674636323047e0, sa,ct);
run_gsw_test(:sigma1, 29.434338510752923e0, sa,ct);
run_gsw_test(:sigma2, 33.609842926904093e0, sa,ct);
run_gsw_test(:sigma3, 37.695147569371784e0, sa,ct);
run_gsw_test(:sigma4, 41.693064726656303e0, sa,ct);
run_gsw_test(:sound_speed, 1527.2011773569989e0, sa,ct,p);
run_gsw_test(:kappa, 4.177024873349404e-010, sa,ct,p);
run_gsw_test(:cabbeling, 9.463053321129075e-6, sa,ct,p);
run_gsw_test(:thermobaric, 1.739078662082863e-12, sa,ct,p);
run_gsw_test(:internal_energy, 79740.482561720783e0, sa,ct,p);
run_gsw_test(:enthalpy, 82761.872939932495e0, sa,ct,p);
run_gsw_test(:dynamic_enthalpy, 2924.5137975399025e0, sa,ct,p);
ρ = rho(sa,ct,p);
run_gsw_test(:sa_from_rho, sa, ρ,ct,p);

# Water column properties, based on the 48-term expression for density
let
    nz = 3
    sa_profile = [35.5e0, 35.7e0, 35.6e0]
    ct_profile = [12.5e0, 15e0, 10e0]
    p_profile = [00e0, 50e0, 100e0]
    lat_profile = [10e0, 10e0, 10e0]

    n2 = zeros(Float64, 2)
    n2_error = zeros(Float64, 2)
    p_mid_n2 = zeros(Float64, 2)
    p_mid_n2_error = zeros(Float64, 2)

    nsquared(sa_profile, ct_profile, p_profile, lat_profile, nz, n2, p_mid_n2);
    n2_error[1] = abs(n2[1] + 0.070960392693051e-3)
    n2_error[2] = abs(n2[2] - 0.175435821615983e-3)
    p_mid_n2_error[1] = abs(p_mid_n2[1] - 25e0);
    p_mid_n2_error[2] = abs(p_mid_n2[2] - 75e0);
    @test n2_error[1] < testTolerances[:n2]
    @test n2_error[1] < testTolerances[:n2]
    @test p_mid_n2_error[1] < testTolerances[:p_mid_n2]
    @test p_mid_n2_error[2] < testTolerances[:p_mid_n2]
end

let
    sa_profile = [35.5e0, 35.7e0, 35.6e0]
    ct_profile = [12.5e0, 15e0, 10e0]
    p_profile = [00e0, 50e0, 100e0]
    nz = 3

    tu = zeros(Float64, 2)
    tu_error = zeros(Float64, 2)
    rsubrho = zeros(Float64, 2)
    rsubrho_error = zeros(Float64, 2)
    p_mid_tursr = zeros(Float64, 2)
    p_mid_tursr_error = zeros(Float64, 2)

    turner_rsubrho(sa_profile, ct_profile, p_profile, nz, tu, rsubrho, p_mid_tursr);
    tu_error[1] = abs(tu[1] + 1.187243981606485e2);
    tu_error[2] = abs(tu[2] - 0.494158257088517e2);
    rsubrho_error[1] = abs(rsubrho[1] - 3.425146897090065e0);
    rsubrho_error[2] = abs(rsubrho[2] - 12.949399443139164e0);
    p_mid_tursr_error[1] = abs(p_mid_tursr[1] - 25e0);
    p_mid_tursr_error[2] = abs(p_mid_tursr[2] - 75e0);
    @test tu_error[1] < testTolerances[:tu]
    @test tu_error[2] < testTolerances[:tu]
    @test rsubrho_error[1] < testTolerances[:rsubrho]
    @test rsubrho_error[2] < testTolerances[:rsubrho]
    @test p_mid_tursr_error[1] < testTolerances[:p_mid_tursr]
    @test p_mid_tursr_error[2] < testTolerances[:p_mid_tursr]
end

let
    sa_profile = [35.5e0, 35.7e0, 35.6e0]
    ct_profile = [12.5e0, 15e0, 10e0]
    p_profile = [00e0, 50e0, 100e0]
    nz = 3

    ipvfn2 = zeros(Float64, 2)
    ipvfn2_error = zeros(Float64, 2)
    p_mid_ipvfn2 = zeros(Float64, 2)
    p_mid_ipvfn2_error = zeros(Float64, 2)

    ipv_vs_fnsquared_ratio(sa_profile, ct_profile, p_profile, nz, ipvfn2, p_mid_ipvfn2)
    ipvfn2_error[1] = abs(ipvfn2[1] - 0.996783975249010e0);
    ipvfn2_error[2] = abs(ipvfn2[2] - 0.992112251478320e0);
    p_mid_ipvfn2_error[1] = abs(p_mid_ipvfn2[1] - 25e0);
    p_mid_ipvfn2_error[2] = abs(p_mid_ipvfn2[2] - 75e0);
    @test ipvfn2_error[1] < testTolerances[:ipvfn2]
    @test ipvfn2_error[2] < testTolerances[:ipvfn2]
    @test p_mid_ipvfn2_error[1] < testTolerances[:p_mid_ipvfn2]
    @test p_mid_ipvfn2_error[2] < testTolerances[:p_mid_ipvfn2]
end

# Freezing temperatures
run_gsw_test(:ct_freezing, -2.1801450326174852e0, sa,p,saturation_fraction);
run_gsw_test(:t_freezing, -2.1765521998023516e0, sa,p,saturation_fraction);

# Isobaric melting enthalpy and isobaric evaporation enthalpy
run_gsw_test(:latentheat_melting, 329330.54839618353e0, sa,p);
run_gsw_test(:latentheat_evap_ct, 2450871.0228523901e0, sa,ct);
run_gsw_test(:latentheat_evap_t, 2462848.2895522709e0, sa,t);

# Basic thermodynamic properties in terms of in-situ t, based on the exact Gibbs function
run_gsw_test(:rho_t_exact, 1027.7128170207150e0, sa,t,p);
run_gsw_test(:pot_rho_t_exact, 1026.8362655887486e0, sa,t,p,p_ref);
run_gsw_test(:alpha_wrt_t_exact, 2.19066952410728916e-4, sa,t,p);
run_gsw_test(:beta_const_t_exact, 7.44744841648729426e-4, sa,t,p);
run_gsw_test(:specvol_t_exact, 9.73034473676164815e-4, sa,t,p);
run_gsw_test(:sound_speed_t_exact, 1512.2053940303056e0, sa,t,p);
run_gsw_test(:kappa_t_exact, 4.25506953386609075e-010, sa,t,p);
run_gsw_test(:enthalpy_t_exact, 62520.680485510929e0, sa,t,p);
run_gsw_test(:cp_t_exact, 3982.7832563441461e0, sa,t,p);

# Library functions of the GSW toolbox
run_gsw_test(:deltasa_atlas, 3.87660373016291727e-3, p,lon,lat);
run_gsw_test(:fdelta, 1.49916256924158942e-004, p,lon,lat);
run_gsw_test(:sa_from_sp_baltic, 35.666154857142850e0, sp,long_bs,lat_bs);
run_gsw_test(:sp_from_sa_baltic, 35.533769845749660e0, sa,long_bs,lat_bs);

