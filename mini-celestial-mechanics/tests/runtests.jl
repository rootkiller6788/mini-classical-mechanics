#!/usr/bin/env julia
include("../src/CelestialMechanics.jl"); using .CelestialMechanics
p=0;f=0
macro t(n,ex) quote try if $(esc(ex)) global p+=1;println("  ✓ ",$(esc(n)))
else global f+=1;println("  ✗ ",$(esc(n)));end catch e global f+=1;println("  ✗ ",$(esc(n)),": ",e);end end end
function main()
println("="^60); println("  Celestial Mechanics Tests"); println("="^60)

# --- Kepler equation ---
println("\n--- Kepler Equation ---")
E=solve_kepler(1.0,0.5); @t("Kepler M=E-e*sinE", abs(E-0.5*sin(E)-1.0)<1e-10)
@t("Kepler E in [0,2pi]", 0<=E<=2pi)
nu=mean_to_true_anomaly(1.0,0.5); @t("Mean→True anomaly", nu isa Float64)
E2=true_to_eccentric_anomaly(nu,0.5); @t("True→Eccentric anomaly", abs(E-E2)<1e-10)
M_round=eccentric_to_mean_anomaly(E2,0.5); @t("Roundtrip M", abs(M_round-1.0)<1e-10)
nu2=eccentric_to_true_anomaly(E,0.5); @t("Roundtrip nu", abs(nu-nu2)<1e-10)

# --- Elements ↔ State ---
println("\n--- Elements ↔ State ---")
a,e,i,Om,om,nu_val=1.0,0.1,0.5,1.0,0.3,1.2
r,v=orbital_elements_to_state(a,e,i,Om,om,nu_val)
@t("r is 3-vector", length(r)==3)
@t("v is 3-vector", length(v)==3)
el2=state_to_orbital_elements(r,v)
@t("Roundtrip a", abs(el2.a-a)<1e-8)
@t("Roundtrip e", abs(el2.e-e)<1e-8)
err=kepler_roundtrip_error(OrbitalElements(a,e,i,Om,om,nu_val))
@t("Roundtrip error small", maximum(values(err))<1e-6)

# --- Two-body ---
println("\n--- Two-body ---")
T=orbital_period(1.0); @t("Period ~2pi", abs(T-2pi)<1e-8)
n=mean_motion(1.0); @t("n*T=2pi", abs(n*T-2pi)<1e-10)
v_esc=escape_velocity(1.0); @t("v_esc=√2", abs(v_esc-sqrt(2))<1e-10)
v_circ=circular_velocity(1.0); @t("v_circ=1", abs(v_circ-1.0)<1e-10)
dt=time_of_flight_elliptic(1.0,0.5,0.0,pi)
@t("TOF > 0", dt>0)
r_p=periapsis_distance(1.0,0.5); @t("r_p=0.5", abs(r_p-0.5)<1e-10)
r_a=apoapsis_distance(1.0,0.5); @t("r_a=1.5", abs(r_a-1.5)<1e-10)
nu_prop=propagate_anomaly(1.0,0.0,0.0,pi)
@t("Propagate half-orb: nu~pi", abs(nu_prop-pi)<1e-8)

# --- Types ---
println("\n--- Types ---")
el=OrbitalElements(1.0,0.1,0.5,1.0,0.3,1.2)
@t("OrbitalElements creation", el.a==1.0)
@t("Elliptic type", orbit_type(el)==:elliptic)
el_hyper=OrbitalElements(-2.0,1.5,0.0,0.0,0.0,0.0)
@t("Hyperbolic type", orbit_type(el_hyper)==:hyperbolic)
ko=KeplerOrbit(el)
@t("KeplerOrbit T>0", ko.T > 0)
@t("KeplerOrbit r_p<r_a", ko.r_p < ko.r_a)

# --- Perturbations ---
println("\n--- Perturbations ---")
a_j2=j2_acceleration([7000.0,0.0,0.0])
@t("J2 accel 3-vector", length(a_j2)==3)
@t("J2 accel non-zero", norm(a_j2)>1e-12)
rates=j2_secular_rates(7000.0,0.001,0.5)
@t("J2 rates Omega_dot<0", rates.Omega_dot < 0)
@t("J2 rates omega_dot exists", rates.omega_dot isa Float64)
a_drag=atmospheric_drag_acceleration([7000.0,0.0,0.0],[0.0,7.5,0.0])
@t("Drag accel 3-vector", length(a_drag)==3)
a_srp=srp_acceleration([7000.0,0.0,0.0],[1.5e8,0.0,0.0])
@t("SRP accel 3-vector", length(a_srp)==3)
a_3b=third_body_acceleration([7000.0,0.0,0.0],[384400.0,0.0,0.0],4902.8)
@t("3rd body accel 3-vector", length(a_3b)==3)
gr=gr_precession_rate(0.387,0.2056)
@t("GR precession >0", gr>0)

# --- Three-body ---
println("\n--- Three-body ---")
pts=lagrange_points(0.01)
@t("5 Lagrange points exist", pts isa LagrangePointData)
CJs=lagrange_jacobi_constants(0.01)
@t("5 Jacobi constants", length(CJs)==5)
hr=hill_region(CJs[1]+0.1,0.01)
@t("Hill region returns string", hr isa String)
soi=sphere_of_influence(1.0,3e-6,1.0); @t("SOI>0", soi>0)
hs=hill_sphere(1.0,3e-6,1.0); @t("Hill>0", hs>0)
rr=roche_limit_rigid(6371.0,5515.0,3340.0); @t("Roche rigid>0", rr>0)
rf=roche_limit_fluid(6371.0,5515.0,3340.0); @t("Roche fluid>rigid", rf>rr)
state0=[0.5,0.0,0.0,0.0,0.5,0.0]
traj_crtbp=integrate_crtbp(state0,0.01,0.1,0.01)
@t("CRTBP integration works", length(traj_crtbp)==11)

# --- Mission ---
println("\n--- Mission ---")
hoh=hohmann_transfer(1.0,1.524)
@t("Hohmann dv>0", hoh.delta_v_total>0)
@t("Hohmann has 2 burns", hoh.delta_v1>0 && hoh.delta_v2>0)
@t("Hohmann a_trans", hoh.a_trans==(1.0+1.524)/2.0)
bi=bi_elliptic_transfer(1.0,1.524,2.0)
@t("Bi-elliptic has 3 burns", bi.delta_v1>0 && bi.delta_v2>0 && bi.delta_v3>0)
dv_plane=simple_plane_change_delta_v(1.0,0.5)
@t("Plane change dv>0", dv_plane>0)
delta=gravity_assist_turn_angle(3.0,6378.0)
@t("Gravity assist turn angle >0", delta>0)
phi=phasing_angle_hohmann(1.0,1.524)
@t("Phasing angle in [0,2pi]", 0<=phi<=2pi)
syn=synodic_period(365.25,687.0)
@t("Synodic period>0", syn>0)

tot=p+f; println("\n"*"="^60); println("  $p/$tot passed"); println(f==0?"ALL PASSED":"$f FAILED")
end; main(); exit(f==0?0:1)
