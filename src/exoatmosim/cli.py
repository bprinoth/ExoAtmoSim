"""Command-line stages for ExoAtmoSim. See `exoatmosim <stage> --help`."""
from __future__ import annotations
import argparse
from pathlib import Path


def main(argv=None):
    parser = argparse.ArgumentParser(prog='exoatmosim',
                        description='Atmosphere, transmission spectrum and instrument observation simulation')
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('simulate', add_help=False, help='Run the complete NASA → FastChem → POSEIDON → instrument workflow')
    config = sub.add_parser('config', help='Create or reuse a NASA planet configuration')
    config.add_argument('--planet', required=True)
    config.add_argument('--instrument', default='CRIRES_K')
    config.add_argument('--data-dir', type=Path, default=Path('data'))
    config.add_argument('--site-config', type=Path)
    config.add_argument('--refresh-config', action='store_true')
    config.add_argument('--fresh', action='store_true')
    atmosphere = sub.add_parser('atmosphere', help='Compute atmospheric T-P and abundance diagnostics')
    atmosphere.add_argument('config', type=Path)
    atmosphere.add_argument('--seed', type=int, default=42)
    atmosphere.add_argument('--overwrite', action='store_true')
    template = sub.add_parser('template', help='Compute a single-species cross-correlation transmission template')
    template.add_argument('config', type=Path)
    template.add_argument('--species', required=True, help='Species name as used in atmospheric abundances (e.g. O3)')
    template.add_argument('--seed', type=int, default=42)
    template.add_argument('--overwrite', action='store_true')
    spectrum = sub.add_parser('spectrum', help='Compute a planetary transmission spectrum')
    spectrum.add_argument('config', type=Path)
    spectrum.add_argument('--backend', default='poseidon', choices=['poseidon'])
    spectrum.add_argument('--seed', type=int, default=42)
    spectrum.add_argument('--overwrite', action='store_true')
    cc = sub.add_parser('cc', help='Prepare and cross-correlate existing observed spectra')
    cc.add_argument('config', type=Path)
    cc.add_argument('--seed', type=int, default=42)
    cc.add_argument('--sysrem', type=int, default=5)
    cc.add_argument('--species', nargs='+', required=True,
                    help='Existing species templates, e.g. O2 O3 H2O')
    cc_mode = cc.add_mutually_exclusive_group()
    cc_mode.add_argument('--prepare-only', action='store_true')
    cc_mode.add_argument('--cc-only', action='store_true')
    cc.add_argument('--overwrite', action='store_true')
    observe = sub.add_parser('observe', help='Observe an existing planetary FITS with constant-S/N orders')
    observe.add_argument('config', type=Path)
    observe.add_argument('--seed', type=int, default=42)
    observe.add_argument('--snr', type=float, default=None,
                         help='Constant S/N for non-VIOLA; VIOLA uses ETC-derived electron noise')
    observe.add_argument('--stellar-npz', type=Path)
    observe.add_argument('--overwrite', action='store_true')
    args, remaining = parser.parse_known_args(argv)
    if args.command == 'simulate':
        from .full_crires import main as full_main
        return full_main(remaining)
    if remaining:
        parser.error('unrecognized arguments: ' + ' '.join(remaining))
    from . import workflow
    if args.command == 'config':
        result = workflow.make_config(args.planet, instrument=args.instrument,
                data_dir=args.data_dir, site_config=args.site_config,
                refresh=args.refresh_config, fresh=args.fresh)
    elif args.command == 'atmosphere':
        result = workflow.calculate_atmosphere(args.config, seed=args.seed, overwrite=args.overwrite)
    elif args.command == 'template':
        result = workflow.calculate_template(args.config, species=args.species,
                                             seed=args.seed, overwrite=args.overwrite)
    elif args.command == 'cc':
        from .cc import cross_correlate_run
        result = cross_correlate_run(args.config, seed=args.seed,
                sysrem=args.sysrem, species=tuple(args.species),
                prepare_only=args.prepare_only, cc_only=args.cc_only,
                overwrite=args.overwrite)
    elif args.command == 'spectrum':
        result = workflow.calculate_spectrum(args.config, backend=args.backend,
                                             seed=args.seed, overwrite=args.overwrite)
    else:
        result = workflow.simulate_observation(args.config, seed=args.seed,
                snr=args.snr, stellar_npz=args.stellar_npz, overwrite=args.overwrite)
    print(f'[SUCCESS] {args.command}: {result}')
    return 0


if __name__ == '__main__':
    main()
