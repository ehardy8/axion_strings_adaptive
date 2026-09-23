(* ::Package:: *)

(* Mathematica port of this session's core-breathing/handoff-smoothing
   diagnostics (2026-09-24, with the user: "it'll be useful to have
   mathematica as well" alongside the checked-in Python pipeline,
   axion_analysis.py/plot_run_summary.py). Deliberately lighter-weight
   than the Python version -- no Moore-switch-aware log(m_r/H) axis, no
   f_a-normalised physical-unit spectrum -- see plot_run_summary.py for
   that. This covers exactly the plots this session's investigation
   actually used: xi(tau), the radial (Higgs) energy components that are
   the core-breathing signature, and the axion spectrum vs k/H at several
   snapshots.

   Usage: set networkFile/spectrumFile below (or pass as script args via
   $ScriptCommandLine), then evaluate the whole file. Column orders match
   AxionStringsLevel.cpp's write_header_line calls exactly -- if that
   header layout ever changes again, update the two association literals
   below (same maintenance note as axion_analysis.py's _NETWORK_COLUMNS/
   _SPECTRUM_COLUMNS). *)

networkFile = "network_scalars.dat";
spectrumFile = "axion_spectrum.dat";

(* -------------------------------------------------------------------- *)
(* network_scalars.dat *)
(* -------------------------------------------------------------------- *)
net = Import[networkFile, "Table", "HeaderLines" -> 1];

netCol = <|"tau" -> 1, "N_p" -> 2, "N_p_weighted" -> 3, "xi" -> 4,
   "xi_weighted" -> 5, "m_r_over_H" -> 6, "rho_tot_unscreened" -> 7,
   "rho_tot_screened" -> 8, "rho_axion_kin_unscreened" -> 9,
   "rho_axion_kin_screened" -> 10, "rho_axion_grad_unscreened" -> 11,
   "rho_axion_grad_screened" -> 12, "rho_radial_kin_unscreened" -> 13,
   "rho_radial_kin_screened" -> 14, "rho_radial_grad_unscreened" -> 15,
   "rho_radial_grad_screened" -> 16, "rho_radial_mass_unscreened" -> 17,
   "rho_radial_mass_screened" -> 18, "n_total" -> 19, "n_unmasked" -> 20,
   "mean_gamma_sq_v_sq" -> 21, "mean_gamma" -> 22,
   "n_velocity_corners" -> 23, "tension_core_only" -> 24,
   "tension_core_plus_tail" -> 25|>;
nc[name_] := net[[All, netCol[name]]];

tau = nc["tau"];

xiPlot = ListLinePlot[Transpose[{tau, nc["xi"]}],
   PlotMarkers -> Automatic, Frame -> True,
   FrameLabel -> {"\[Tau]", "\[Xi]"}, PlotLabel -> "String network density",
   GridLines -> Automatic, ImageSize -> 480];

radialPlot = ListLinePlot[
   {Transpose[{tau, nc["rho_radial_kin_unscreened"]}],
    Transpose[{tau, nc["rho_radial_grad_unscreened"]}],
    Transpose[{tau, nc["rho_radial_mass_unscreened"]}]},
   PlotLegends -> {"kinetic", "gradient", "mass"},
   PlotMarkers -> Automatic, Frame -> True,
   FrameLabel -> {"\[Tau]", "\[Rho] (unscreened)"},
   PlotLabel -> "Radial (Higgs) energy -- core-breathing signature",
   GridLines -> Automatic, ImageSize -> 480];

npPlot = ListLinePlot[Transpose[{tau, nc["N_p"]}],
   PlotMarkers -> Automatic, Frame -> True,
   FrameLabel -> {"\[Tau]", "N_p"}, PlotLabel -> "Plaquette count",
   GridLines -> Automatic, ImageSize -> 480];

(* -------------------------------------------------------------------- *)
(* axion_spectrum.dat *)
(* -------------------------------------------------------------------- *)
spec = Import[spectrumFile, "Table", "HeaderLines" -> 1];

specCol = <|"tau" -> 1, "mode_index" -> 2, "k_over_H" -> 3,
   "shell_average_screened" -> 4, "shell_average_unscreened" -> 5|>;

specTaus = Union[spec[[All, specCol["tau"]]]];
nPick = Min[6, Length[specTaus]];
pickIdx = If[nPick <= 1, {1}, Round[Subdivide[1, Length[specTaus], nPick - 1]]];
picks = specTaus[[pickIdx]];

spectrumCurve[t_, which_] := Module[{rows},
   rows = Select[spec, #[[specCol["tau"]]] == t &];
   Select[Transpose[{rows[[All, specCol["k_over_H"]]],
      rows[[All, specCol[which]]]}], #[[1]] > 0 && #[[2]] > 0 &]];

spectrumPlot = ListLogLogPlot[
   spectrumCurve[#, "shell_average_unscreened"] & /@ picks,
   PlotLegends -> ("\[Tau]=" <> ToString[NumberForm[#, 4]] & /@ picks),
   Frame -> True, FrameLabel -> {"k/H", "shell_average (unscreened)"},
   PlotLabel -> "Axion spectrum vs k/H at several times",
   GridLines -> Automatic, ImageSize -> 480];

(* Exported separately, not combined into one GraphicsGrid -- combining
   into a fixed-size grid squeezes/clips the legends on the radial-energy
   and spectrum panels. Matches plot_run_summary.py's own convention of
   one file per plot. *)
Export["xi_vs_tau.png", xiPlot];
Export["N_p_vs_tau.png", npPlot];
Export["radial_energy_vs_tau.png", radialPlot];
Export["axion_spectrum_vs_k_over_H.png", spectrumPlot];
Print["Wrote xi_vs_tau.png, N_p_vs_tau.png, radial_energy_vs_tau.png, ",
  "axion_spectrum_vs_k_over_H.png"];
