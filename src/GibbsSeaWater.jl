module GibbsSeaWater
import Base.beta

export add_barrier, add_mean, adiabatic_lapse_rate_from_ct, alpha,
       alpha_on_beta, alpha_wrt_t_exact, beta_const_t_exact, beta, c_from_sp,
       cabbeling, ct_freezing, ct_from_pt, ct_from_t, deltasa_atlas,
       deltasa_from_sp, dynamic_enthalpy, enthalpy, enthalpy_sso_0_p,
       enthalpy_t_exact, cp_t_exact, entropy_from_t, entropy_part,
       entropy_part_zerop, fdelta, gibbs, gibbs_pt0_pt0, grav,
       hill_ratio_at_sp2, indx, internal_energy, ipv_vs_fnsquared_ratio, kappa,
       kappa_t_exact, latentheat_evap_ct, latentheat_evap_t,
       latentheat_melting, nsquared, pot_rho_t_exact, pt0_from_t, pt_from_ct,
       pt_from_t, rho, rho_first_derivatives, rho_t_exact, saar, sa_from_rho,
       sa_from_sp_baltic, sa_from_sp, sa_from_sstar, sigma0, sigma1, sigma2,
       sigma3, sigma4, sound_speed, sound_speed_t_exact, specvol_anom, specvol,
       specvol_sso_0_p, specvol_t_exact, sp_from_c, sp_from_sa_baltic,
       sp_from_sa, sp_from_sk, sp_from_sr, sp_from_sstar, sr_from_sp,
       sstar_from_sa, sstar_from_sp, t_freezing, t_from_ct, thermobaric,
       turner_rsubrho, xinterp1, z_from_p

ext = @linux? (".so" : @osx? ( ".dylib" : @windows? ( ".dll" : "" )))
if ext == ""
    error("Platform not linux, OS X, or Windows")
end

const _MODULEPATH = Pkg.dir("GibbsSeaWater")
const _LIBGSWPATH = joinpath(_MODULEPATH, "deps/build/libgswteos-10"*ext)

function add_barrier(input_data, lon, lat, long_grid, lat_grid, dlong_grid, dlat_grid, output_data)
    return ccall((:gsw_add_barrier, _LIBGSWPATH), Void, (Ptr{Float64}, Float64, Float64, Float64, Float64, Float64, Float64, Ptr{Float64}), pointer(input_data), lon, lat, long_grid, lat_grid, dlong_grid, dlat_grid, pointer(output_data))
end

function add_mean(data_in, lon, lat, data_out)
    return ccall((:gsw_add_mean, _LIBGSWPATH), Void, (Ptr{Float64}, Float64, Float64, Ptr{Float64}), pointer(data_in), lon, lat, pointer(data_out))
end

function adiabatic_lapse_rate_from_ct(sa, ct, p)
    return ccall((:gsw_adiabatic_lapse_rate_from_ct, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, ct, p)
end

function alpha(sa, ct, p)
    return ccall((:gsw_alpha, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, ct, p)
end

function alpha_on_beta(sa, ct, p)
    return ccall((:gsw_alpha_on_beta, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, ct, p)
end

function alpha_wrt_t_exact(sa, t, p)
    return ccall((:gsw_alpha_wrt_t_exact, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, t, p)
end

function beta_const_t_exact(sa, t, p)
    return ccall((:gsw_beta_const_t_exact, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, t, p)
end

function beta(sa, ct, p)
    return ccall((:gsw_beta, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, ct, p)
end

function c_from_sp(sp, t, p)
    return ccall((:gsw_c_from_sp, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sp, t, p)
end

function cabbeling(sa, ct, p)
    return ccall((:gsw_cabbeling, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, ct, p)
end

function ct_freezing(sa, p, saturation_fraction)
    return ccall((:gsw_ct_freezing, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, p, saturation_fraction)
end

function ct_from_pt(sa, pt)
    return ccall((:gsw_ct_from_pt, _LIBGSWPATH), Float64, (Float64, Float64), sa, pt)
end

function ct_from_t(sa, t, p)
    return ccall((:gsw_ct_from_t, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, t, p)
end

function deltasa_atlas(p, lon, lat)
    return ccall((:gsw_deltasa_atlas, _LIBGSWPATH), Float64, (Float64, Float64, Float64), p, lon, lat)
end

function deltasa_from_sp(sp, p, lon, lat)
    return ccall((:gsw_deltasa_from_sp, _LIBGSWPATH), Float64, (Float64, Float64, Float64, Float64), sp, p, lon, lat)
end

function dynamic_enthalpy(sa, ct, p)
    return ccall((:gsw_dynamic_enthalpy, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, ct, p)
end

function enthalpy(sa, ct, p)
    return ccall((:gsw_enthalpy, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, ct, p)
end

function enthalpy_sso_0_p(p)
    return ccall((:gsw_enthalpy_sso_0_p, _LIBGSWPATH), Float64, (Float64,), p)
end

function enthalpy_t_exact(sa, t, p)
    return ccall((:gsw_enthalpy_t_exact, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, t, p)
end

function cp_t_exact(sa, t, p)
    return ccall((:gsw_cp_t_exact, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, t, p)
end

function entropy_from_t(sa, t, p)
    return ccall((:gsw_entropy_from_t, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, t, p)
end

function entropy_part(sa, t, p)
    return ccall((:gsw_entropy_part, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, t, p)
end

function entropy_part_zerop(sa, pt0)
    return ccall((:gsw_entropy_part_zerop, _LIBGSWPATH), Float64, (Float64, Float64), sa, pt0)
end

function fdelta(p, lon, lat)
    return ccall((:gsw_fdelta, _LIBGSWPATH), Float64, (Float64, Float64, Float64), p, lon, lat)
end

function gibbs(ns, nt, np, sa, t, p)
    return ccall((:gsw_gibbs, _LIBGSWPATH), Float64, (Int64, Int64, Int64, Float64, Float64, Float64), ns, nt, np, sa, t, p)
end

function gibbs_pt0_pt0(sa, pt0)
    return ccall((:gsw_gibbs_pt0_pt0, _LIBGSWPATH), Float64, (Float64, Float64), sa, pt0)
end

function grav(lat, p)
    return ccall((:gsw_grav, _LIBGSWPATH), Float64, (Float64, Float64), lat, p)
end

function hill_ratio_at_sp2(t)
    return ccall((:gsw_hill_ratio_at_sp2, _LIBGSWPATH), Float64, (Float64,), t)
end

function indx(x, n, z)
    return ccall((:gsw_indx, _LIBGSWPATH), Int64, (Ptr{Float64}, Int64, Float64), pointer(x), n, z)
end

function internal_energy(sa, ct, p)
    return ccall((:gsw_internal_energy, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, ct, p)
end

function ipv_vs_fnsquared_ratio(sa, ct, p, nz, ipv_vs_fnsquared_ratio, p_mid)
    return ccall((:gsw_ipv_vs_fnsquared_ratio, _LIBGSWPATH), Void, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64, Ptr{Float64}, Ptr{Float64}), pointer(sa), pointer(ct), pointer(p), nz, pointer(ipv_vs_fnsquared_ratio), pointer(p_mid))
end

function kappa(sa, ct, p)
    return ccall((:gsw_kappa, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, ct, p)
end

function kappa_t_exact(sa, t, p)
    return ccall((:gsw_kappa_t_exact, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, t, p)
end

function latentheat_evap_ct(sa, ct)
    return ccall((:gsw_latentheat_evap_ct, _LIBGSWPATH), Float64, (Float64, Float64), sa, ct)
end

function latentheat_evap_t(sa, t)
    return ccall((:gsw_latentheat_evap_t, _LIBGSWPATH), Float64, (Float64, Float64), sa, t)
end

function latentheat_melting(sa, p)
    return ccall((:gsw_latentheat_melting, _LIBGSWPATH), Float64, (Float64, Float64), sa, p)
end

function nsquared(sa, ct, p, lat, nz, n2, p_mid)
    return ccall((:gsw_nsquared, _LIBGSWPATH), Void, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64, Ptr{Float64}, Ptr{Float64}), pointer(sa), pointer(ct), pointer(p), pointer(lat), nz, pointer(n2), pointer(p_mid))
end

function pot_rho_t_exact(sa, t, p, p_ref)
    return ccall((:gsw_pot_rho_t_exact, _LIBGSWPATH), Float64, (Float64, Float64, Float64, Float64), sa, t, p, p_ref)
end

function pt0_from_t(sa, t, p)
    return ccall((:gsw_pt0_from_t, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, t, p)
end

function pt_from_ct(sa, ct)
    return ccall((:gsw_pt_from_ct, _LIBGSWPATH), Float64, (Float64, Float64), sa, ct)
end

function pt_from_t(sa, t, p, p_ref)
    return ccall((:gsw_pt_from_t, _LIBGSWPATH), Float64, (Float64, Float64, Float64, Float64), sa, t, p, p_ref)
end

function rho(sa, ct, p)
    return ccall((:gsw_rho, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, ct, p)
end

function rho_first_derivatives(sa, ct, p, drho_dsa, drho_dct, drho_dp)
    return ccall((:gsw_rho_first_derivatives, _LIBGSWPATH), Void, (Float64, Float64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), sa, ct, p, pointer(drho_dsa), pointer(drho_dct), pointer(drho_dp))
end

function rho_t_exact(sa, t, p)
    return ccall((:gsw_rho_t_exact, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, t, p)
end

function saar(p, lon, lat)
    return ccall((:gsw_saar, _LIBGSWPATH), Float64, (Float64, Float64, Float64), p, lon, lat)
end

function sa_from_rho(rho, ct, p)
    return ccall((:gsw_sa_from_rho, _LIBGSWPATH), Float64, (Float64, Float64, Float64), rho, ct, p)
end

function sa_from_sp_baltic(sp, lon, lat)
    return ccall((:gsw_sa_from_sp_baltic, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sp, lon, lat)
end

function sa_from_sp(sp, p, lon, lat)
    return ccall((:gsw_sa_from_sp, _LIBGSWPATH), Float64, (Float64, Float64, Float64, Float64), sp, p, lon, lat)
end

function sa_from_sstar(sstar, p, lon, lat)
    return ccall((:gsw_sa_from_sstar, _LIBGSWPATH), Float64, (Float64, Float64, Float64, Float64), sstar, p, lon, lat)
end

function sigma0(sa, ct)
    return ccall((:gsw_sigma0, _LIBGSWPATH), Float64, (Float64, Float64), sa, ct)
end

function sigma1(sa, ct)
    return ccall((:gsw_sigma1, _LIBGSWPATH), Float64, (Float64, Float64), sa, ct)
end

function sigma2(sa, ct)
    return ccall((:gsw_sigma2, _LIBGSWPATH), Float64, (Float64, Float64), sa, ct)
end

function sigma3(sa, ct)
    return ccall((:gsw_sigma3, _LIBGSWPATH), Float64, (Float64, Float64), sa, ct)
end

function sigma4(sa, ct)
    return ccall((:gsw_sigma4, _LIBGSWPATH), Float64, (Float64, Float64), sa, ct)
end

function sound_speed(sa, ct, p)
    return ccall((:gsw_sound_speed, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, ct, p)
end

function sound_speed_t_exact(sa, t, p)
    return ccall((:gsw_sound_speed_t_exact, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, t, p)
end

function specvol_anom(sa, ct, p)
    return ccall((:gsw_specvol_anom, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, ct, p)
end

function specvol(sa, ct, p)
    return ccall((:gsw_specvol, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, ct, p)
end

function specvol_sso_0_p(p)
    return ccall((:gsw_specvol_sso_0_p, _LIBGSWPATH), Float64, (Float64,), p)
end

function specvol_t_exact(sa, t, p)
    return ccall((:gsw_specvol_t_exact, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, t, p)
end

function sp_from_c(c, t, p)
    return ccall((:gsw_sp_from_c, _LIBGSWPATH), Float64, (Float64, Float64, Float64), c, t, p)
end

function sp_from_sa_baltic(sa, lon, lat)
    return ccall((:gsw_sp_from_sa_baltic, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, lon, lat)
end

function sp_from_sa(sa, p, lon, lat)
    return ccall((:gsw_sp_from_sa, _LIBGSWPATH), Float64, (Float64, Float64, Float64, Float64), sa, p, lon, lat)
end

function sp_from_sk(sk)
    return ccall((:gsw_sp_from_sk, _LIBGSWPATH), Float64, (Float64,), sk)
end

function sp_from_sr(sr)
    return ccall((:gsw_sp_from_sr, _LIBGSWPATH), Float64, (Float64,), sr)
end

function sp_from_sstar(sstar, p, lon, lat)
    return ccall((:gsw_sp_from_sstar, _LIBGSWPATH), Float64, (Float64, Float64, Float64, Float64), sstar, p, lon, lat)
end

function sr_from_sp(sp)
    return ccall((:gsw_sr_from_sp, _LIBGSWPATH), Float64, (Float64,), sp)
end

function sstar_from_sa(sa, p, lon, lat)
    return ccall((:gsw_sstar_from_sa, _LIBGSWPATH), Float64, (Float64, Float64, Float64, Float64), sa, p, lon, lat)
end

function sstar_from_sp(sp, p, lon, lat)
    return ccall((:gsw_sstar_from_sp, _LIBGSWPATH), Float64, (Float64, Float64, Float64, Float64), sp, p, lon, lat)
end

function t_freezing(sa, p, saturation_fraction)
    return ccall((:gsw_t_freezing, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, p, saturation_fraction)
end

function t_from_ct(sa, ct, p)
    return ccall((:gsw_t_from_ct, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, ct, p)
end

function thermobaric(sa, ct, p)
    return ccall((:gsw_thermobaric, _LIBGSWPATH), Float64, (Float64, Float64, Float64), sa, ct, p)
end

function turner_rsubrho(sa, ct, p, nz, tu, rsubrho, p_mid)
    return ccall((:gsw_turner_rsubrho, _LIBGSWPATH), Void, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), pointer(sa), pointer(ct), pointer(p), nz, pointer(tu), pointer(rsubrho), pointer(p_mid))
end

function xinterp1(x, y, n, x0)
    return ccall((:gsw_xinterp1, _LIBGSWPATH), Float64, (Ptr{Float64}, Ptr{Float64}, Int64, Float64), pointer(x), pointer(y), n, x0)
end

function z_from_p(p, lat)
    return ccall((:gsw_z_from_p, _LIBGSWPATH), Float64, (Float64, Float64), p, lat)
end


end # module
