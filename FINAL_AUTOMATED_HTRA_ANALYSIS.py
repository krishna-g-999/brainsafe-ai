#!/usr/bin/env python3
"""
================================================================================
FULLY AUTOMATED HtrA STRUCTURE ANALYSIS PIPELINE
================================================================================

Complete structural comparison of HtrA protease catalytic triads across:
- Human: HtrA1 (3TJN), HtrA2 (1LCY), HtrA3 (4RI0)
- Bacterial: H. pylori (5Y2D, 7XS0)

USAGE:
    python3 FINAL_AUTOMATED_HTRA_ANALYSIS.py

OUTPUT:
    - catalytic_triad_distances.csv (CORE DATA)
    - MASTER_ANALYSIS_REPORT.csv
    - ANALYSIS_REPORT.md
    - analysis_results.json
    - Publication-ready figures (if PyMOL available)

RUNTIME: ~4-5 minutes

NO MANUAL STEPS REQUIRED - FULLY AUTOMATED
================================================================================
"""

import os
import sys
import json
import urllib.request
from datetime import datetime
from pathlib import Path

print("="*80)
print("AUTOMATED HtrA STRUCTURE ANALYSIS PIPELINE")
print("="*80)
print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# ============================================================================
# CONFIGURATION
# ============================================================================

class Config:
    """Centralized configuration for all analysis parameters"""
    
    # PDB structure definitions with catalytic triad residue numbers
    PDB_IDS = {
        '3TJN': {
            'name': 'HtrA1_apo',
            'his': 220,
            'asp': 250,
            'ser': 328,
            'chain': 'A',
            'description': 'Human HtrA1 (apo form)'
        },
        '1LCY': {
            'name': 'HtrA2_pyramid',
            'his': 198,
            'asp': 228,
            'ser': 306,
            'chain': 'A',
            'description': 'Human HtrA2 (pyramidal intermediate)'
        },
        '4RI0': {
            'name': 'HtrA3_locked',
            'his': 191,
            'asp': 227,
            'ser': 305,
            'chain': 'A',
            'description': 'Human HtrA3 (allosterically locked)'
        },
        '5Y2D': {
            'name': 'Hpylori_trimer',
            'his': 116,
            'asp': 147,
            'ser': 221,
            'chain': 'A',
            'description': 'H. pylori HtrA (active trimer)'
        },
        '7XS0': {
            'name': 'Hpylori_monomer',
            'his': 116,
            'asp': 147,
            'ser': 221,
            'chain': 'A',
            'description': 'H. pylori HtrA (inactive monomer)'
        }
    }
    
    # Directory structure
    BASE_DIR = 'htra_analysis_automated'
    PDB_DIR = os.path.join(BASE_DIR, 'pdb_files')
    RESULTS_DIR = os.path.join(BASE_DIR, 'results')
    DISTANCE_DIR = os.path.join(RESULTS_DIR, 'distances')
    ALIGNMENT_DIR = os.path.join(RESULTS_DIR, 'alignment')
    CASTP_DIR = os.path.join(PDB_DIR, 'for_castp')
    FIGURES_DIR = os.path.join(BASE_DIR, 'pymol_figures')
    
    # PDB download URL
    PDB_URL = "https://files.rcsb.org/download/{}.pdb"

# ============================================================================
# SETUP
# ============================================================================

def setup_directories():
    """Create all required directories"""
    print("Setting up directory structure...")
    dirs = [
        Config.BASE_DIR,
        Config.PDB_DIR,
        Config.RESULTS_DIR,
        Config.DISTANCE_DIR,
        Config.ALIGNMENT_DIR,
        Config.CASTP_DIR,
        Config.FIGURES_DIR
    ]
    for d in dirs:
        Path(d).mkdir(parents=True, exist_ok=True)
    print(f"✓ Created directories in: {Config.BASE_DIR}\n")

def check_dependencies():
    """Check and install required Python packages"""
    print("Checking dependencies...")
    
    required = {
        'Bio': 'biopython',
        'numpy': 'numpy'
    }
    
    missing = []
    for module, package in required.items():
        try:
            __import__(module)
            print(f"✓ {package} found")
        except ImportError:
            print(f"✗ {package} not found")
            missing.append(package)
    
    if missing:
        print(f"\nInstalling missing packages: {', '.join(missing)}")
        import subprocess
        for package in missing:
            subprocess.check_call([sys.executable, "-m", "pip", "install", package, "-q"])
            print(f"✓ Installed {package}")
    
    print()

# ============================================================================
# PDB FILE HANDLING
# ============================================================================

def download_pdb(pdb_id):
    """Download PDB file from RCSB"""
    output_path = os.path.join(Config.PDB_DIR, f"{pdb_id}.pdb")
    
    if os.path.exists(output_path):
        print(f"  ✓ {pdb_id}.pdb (cached)")
        return output_path
    
    url = Config.PDB_URL.format(pdb_id)
    try:
        print(f"  ↓ Downloading {pdb_id}.pdb...", end='')
        urllib.request.urlretrieve(url, output_path)
        print(" Done")
        return output_path
    except Exception as e:
        print(f" Failed: {e}")
        return None

def download_all_structures():
    """Download all PDB structures"""
    print("Downloading PDB structures...")
    
    pdb_files = {}
    for pdb_id in Config.PDB_IDS.keys():
        path = download_pdb(pdb_id)
        if path:
            pdb_files[pdb_id] = path
    
    print(f"✓ Downloaded {len(pdb_files)}/{len(Config.PDB_IDS)} structures\n")
    return pdb_files

def prepare_for_castp(pdb_id, input_path):
    """Prepare cleaned PDB file for CASTp analysis"""
    output_path = os.path.join(Config.CASTP_DIR, f"{pdb_id}_clean.pdb")
    
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            if line.startswith('ATOM'):
                outfile.write(line)
            elif line.startswith('END'):
                outfile.write(line)
                break
    
    return output_path

# ============================================================================
# STRUCTURAL ANALYSIS
# ============================================================================

def calculate_distance(coord1, coord2):
    """Calculate Euclidean distance between two 3D coordinates"""
    import numpy as np
    return np.sqrt(np.sum((np.array(coord1) - np.array(coord2))**2))

def get_atom_coordinates(pdb_file, chain, residue_num, atom_name):
    """Extract coordinates for specific atom from PDB file"""
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                if (line[21:22].strip() == chain and
                    int(line[22:26].strip()) == residue_num and
                    line[12:16].strip() == atom_name):
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    return [x, y, z]
    return None

def measure_catalytic_triad(pdb_id, pdb_file, config):
    """Measure distances in catalytic triad"""
    chain = config['chain']
    his_num = config['his']
    asp_num = config['asp']
    ser_num = config['ser']
    
    # Get coordinates for catalytic residues
    his_ne2 = get_atom_coordinates(pdb_file, chain, his_num, 'NE2')
    asp_od1 = get_atom_coordinates(pdb_file, chain, asp_num, 'OD1')
    ser_og = get_atom_coordinates(pdb_file, chain, ser_num, 'OG')
    
    if not all([his_ne2, asp_od1, ser_og]):
        print(f"  ✗ Could not find all catalytic triad atoms in {pdb_id}")
        return None
    
    # Calculate distances
    his_asp = calculate_distance(his_ne2, asp_od1)
    his_ser = calculate_distance(his_ne2, ser_og)
    asp_ser = calculate_distance(asp_od1, ser_og)
    
    # Classify catalytic state based on distances
    if his_asp < 3.5 and his_ser < 4.5:
        state = "COMPETENT (Active)"
    elif his_asp > 6.0 or his_ser > 5.5:
        state = "ALLOSTERICALLY LOCKED"
    elif his_asp > 4.5:
        state = "DISTORTED (Inactive)"
    else:
        state = "pH-DEPENDENT"
    
    return {
        'pdb_id': pdb_id,
        'name': config['name'],
        'his_asp': his_asp,
        'his_ser': his_ser,
        'asp_ser': asp_ser,
        'state': state,
        'residues': f"His{his_num}-Asp{asp_num}-Ser{ser_num}"
    }

def analyze_all_structures(pdb_files):
    """Analyze all structures and measure catalytic triads"""
    print("Analyzing catalytic triad distances...")
    
    results = []
    for pdb_id, pdb_file in pdb_files.items():
        config = Config.PDB_IDS[pdb_id]
        print(f"  Measuring {pdb_id} ({config['name']})...")
        
        result = measure_catalytic_triad(pdb_id, pdb_file, config)
        if result:
            results.append(result)
            print(f"    His-Asp: {result['his_asp']:.2f} Å")
            print(f"    His-Ser: {result['his_ser']:.2f} Å")
            print(f"    State: {result['state']}")
    
    print(f"✓ Analyzed {len(results)} structures\n")
    return results

def calculate_rmsd_alignments(pdb_files):
    """Calculate RMSD between all structure pairs"""
    print("Calculating structural alignments (RMSD)...")
    
    try:
        from Bio.PDB import PDBParser, Superimposer
        import numpy as np
    except ImportError:
        print("  ⚠ BioPython not fully available, skipping RMSD calculations\n")
        return []
    
    parser = PDBParser(QUIET=True)
    rmsd_results = []
    
    pdb_ids = list(pdb_files.keys())
    for i, pdb1 in enumerate(pdb_ids):
        for pdb2 in pdb_ids[i+1:]:
            try:
                struct1 = parser.get_structure(pdb1, pdb_files[pdb1])
                struct2 = parser.get_structure(pdb2, pdb_files[pdb2])
                
                chain1 = Config.PDB_IDS[pdb1]['chain']
                chain2 = Config.PDB_IDS[pdb2]['chain']
                
                atoms1 = [atom for atom in struct1[0][chain1].get_atoms() if atom.name == 'CA']
                atoms2 = [atom for atom in struct2[0][chain2].get_atoms() if atom.name == 'CA']
                
                min_len = min(len(atoms1), len(atoms2))
                atoms1 = atoms1[:min_len]
                atoms2 = atoms2[:min_len]
                
                super_imposer = Superimposer()
                super_imposer.set_atoms(atoms1, atoms2)
                rmsd = super_imposer.rms
                
                rmsd_results.append({
                    'structure1': f"{pdb1} ({Config.PDB_IDS[pdb1]['name']})",
                    'structure2': f"{pdb2} ({Config.PDB_IDS[pdb2]['name']})",
                    'rmsd': rmsd,
                    'aligned_residues': min_len
                })
                
                print(f"  {pdb1} vs {pdb2}: RMSD = {rmsd:.2f} Å")
                
            except Exception as e:
                print(f"  ✗ Failed to align {pdb1} vs {pdb2}: {e}")
    
    print(f"✓ Calculated {len(rmsd_results)} pairwise alignments\n")
    return rmsd_results

# ============================================================================
# OUTPUT GENERATION
# ============================================================================

def save_distance_csv(results):
    """Save catalytic triad distances to CSV"""
    output_file = os.path.join(Config.DISTANCE_DIR, 'catalytic_triad_distances.csv')
    
    with open(output_file, 'w') as f:
        f.write("PDB_ID,Protein,Residues,His_Asp_Distance_(Å),His_Ser_Distance_(Å),Asp_Ser_Distance_(Å),Catalytic_State\n")
        for r in results:
            f.write(f"{r['pdb_id']},{r['name']},{r['residues']},{r['his_asp']:.2f},{r['his_ser']:.2f},{r['asp_ser']:.2f},{r['state']}\n")
    
    print(f"✓ Saved: {output_file}")
    return output_file

def save_rmsd_csv(rmsd_results):
    """Save RMSD comparison to CSV"""
    if not rmsd_results:
        return None
    
    output_file = os.path.join(Config.ALIGNMENT_DIR, 'rmsd_comparison.csv')
    
    with open(output_file, 'w') as f:
        f.write("Structure_1,Structure_2,RMSD_(Å),Aligned_Residues\n")
        for r in rmsd_results:
            f.write(f"{r['structure1']},{r['structure2']},{r['rmsd']:.2f},{r['aligned_residues']}\n")
    
    print(f"✓ Saved: {output_file}")
    return output_file

def save_master_report(results):
    """Save comprehensive master report"""
    output_file = os.path.join(Config.BASE_DIR, 'MASTER_ANALYSIS_REPORT.csv')
    
    with open(output_file, 'w') as f:
        f.write("PDB_ID,Protein,Description,Catalytic_Triad,His_Asp_(Å),His_Ser_(Å),Asp_Ser_(Å),State\n")
        for r in results:
            pdb_id = r['pdb_id']
            desc = Config.PDB_IDS[pdb_id]['description']
            f.write(f"{r['pdb_id']},{r['name']},\"{desc}\",{r['residues']},{r['his_asp']:.2f},{r['his_ser']:.2f},{r['asp_ser']:.2f},{r['state']}\n")
    
    print(f"✓ Saved: {output_file}")
    return output_file

def save_json_report(results, rmsd_results):
    """Save analysis results in JSON format"""
    output_file = os.path.join(Config.BASE_DIR, 'analysis_results.json')
    
    data = {
        'timestamp': datetime.now().isoformat(),
        'structures_analyzed': len(results),
        'catalytic_triad_distances': results,
        'rmsd_alignments': rmsd_results,
        'pdb_configurations': Config.PDB_IDS
    }
    
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)
    
    print(f"✓ Saved: {output_file}")
    return output_file

def save_markdown_report(results, rmsd_results):
    """Save comprehensive markdown report"""
    output_file = os.path.join(Config.BASE_DIR, 'ANALYSIS_REPORT.md')
    
    with open(output_file, 'w') as f:
        f.write("# HtrA Catalytic Triad Analysis Report\n\n")
        f.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("---\n\n")
        
        f.write("## Executive Summary\n\n")
        f.write(f"Analyzed {len(results)} HtrA protease structures for catalytic triad geometry.\n\n")
        
        f.write("## Catalytic Triad Distances\n\n")
        f.write("| PDB ID | Protein | His-Asp (Å) | His-Ser (Å) | Asp-Ser (Å) | State |\n")
        f.write("|--------|---------|-------------|-------------|-------------|-------|\n")
        for r in results:
            f.write(f"| {r['pdb_id']} | {r['name']} | {r['his_asp']:.2f} | {r['his_ser']:.2f} | {r['asp_ser']:.2f} | {r['state']} |\n")
        
        f.write("\n## Key Findings\n\n")
        
        active = [r for r in results if 'COMPETENT' in r['state']]
        inactive = [r for r in results if 'DISTORTED' in r['state'] or 'LOCKED' in r['state']]
        
        f.write(f"- **Active/Competent structures:** {len(active)}\n")
        f.write(f"- **Inactive/Distorted structures:** {len(inactive)}\n\n")
        
        f.write("### Catalytic State Classification\n\n")
        f.write("**Competent (Active):**\n")
        for r in active:
            f.write(f"- {r['pdb_id']} ({r['name']}): His-Asp = {r['his_asp']:.2f} Å\n")
        
        f.write("\n**Inactive/Distorted:**\n")
        for r in inactive:
            f.write(f"- {r['pdb_id']} ({r['name']}): His-Asp = {r['his_asp']:.2f} Å\n")
        
        if rmsd_results:
            f.write("\n## Structural Alignments (RMSD)\n\n")
            f.write("| Structure 1 | Structure 2 | RMSD (Å) | Residues |\n")
            f.write("|-------------|-------------|----------|----------|\n")
            for r in rmsd_results:
                f.write(f"| {r['structure1']} | {r['structure2']} | {r['rmsd']:.2f} | {r['aligned_residues']} |\n")
        
        f.write("\n## Next Steps\n\n")
        f.write("1. **CASTp Analysis**: Upload cleaned PDB files from `pdb_files/for_castp/`\n")
        f.write("2. **PDBePISA**: Analyze oligomerization interfaces\n")
        f.write("3. **Publication**: Use figures from `pymol_figures/` directory\n")
        f.write("\n---\n\n")
        f.write("*Analysis completed successfully*\n")
    
    print(f"✓ Saved: {output_file}")
    return output_file

# ============================================================================
# PYMOL FIGURE GENERATION
# ============================================================================

def generate_pymol_figures():
    """Generate publication-quality figures using PyMOL"""
    print("Attempting to generate PyMOL figures...")
    
    try:
        import pymol
        from pymol import cmd
        pymol.finish_launching(['pymol', '-c'])
        
        print("✓ PyMOL available - generating figures...\n")
        
        cmd.set('ray_trace_mode', 1)
        cmd.set('ray_shadow', 'on')
        cmd.set('antialias', 2)
        cmd.bg_color('white')
        
        figures_generated = 0
        
        for pdb_id, config in Config.PDB_IDS.items():
            try:
                print(f"  Rendering {pdb_id} ({config['name']})...")
                
                cmd.fetch(pdb_id, pdb_id)
                cmd.hide('all')
                cmd.show('cartoon', pdb_id)
                cmd.color('white', pdb_id)
                
                triad_sel = f"{pdb_id} and chain {config['chain']} and (resi {config['his']} or resi {config['asp']} or resi {config['ser']})"
                cmd.select('cat_triad', triad_sel)
                cmd.show('sticks', 'cat_triad')
                cmd.color('red', 'cat_triad')
                
                his_atom = f"{pdb_id} and chain {config['chain']} and resi {config['his']} and name NE2"
                asp_atom = f"{pdb_id} and chain {config['chain']} and resi {config['asp']} and name OD1"
                ser_atom = f"{pdb_id} and chain {config['chain']} and resi {config['ser']} and name OG"
                
                cmd.distance('d1', his_atom, asp_atom)
                cmd.distance('d2', his_atom, ser_atom)
                cmd.distance('d3', asp_atom, ser_atom)
                
                cmd.label(f"{pdb_id} and resi {config['his']}", f"'His{config['his']}'")
                cmd.label(f"{pdb_id} and resi {config['asp']}", f"'Asp{config['asp']}'")
                cmd.label(f"{pdb_id} and resi {config['ser']}", f"'Ser{config['ser']}'")
                
                cmd.zoom('cat_triad')
                cmd.orient()
                
                output_file = os.path.join(Config.FIGURES_DIR, f"{pdb_id}_catalytic_triad.png")
                cmd.png(output_file, width=2400, height=2400, dpi=300, ray=1)
                
                print(f"    ✓ Saved: {output_file}")
                figures_generated += 1
                
                cmd.delete('all')
                
            except Exception as e:
                print(f"    ✗ Failed: {e}")
        
        cmd.quit()
        print(f"\n✓ Generated {figures_generated} figures\n")
        return True
        
    except ImportError:
        print("⚠ PyMOL not available - skipping figure generation")
        print("  (Analysis data still complete without figures)\n")
        return False
    except Exception as e:
        print(f"⚠ PyMOL error: {e}")
        print("  (Analysis data still complete without figures)\n")
        return False

# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    """Execute complete automated analysis pipeline"""
    
    print("\n" + "="*80)
    print("STEP 1: SETUP")
    print("="*80 + "\n")
    
    setup_directories()
    check_dependencies()
    
    print("="*80)
    print("STEP 2: DOWNLOAD PDB STRUCTURES")
    print("="*80 + "\n")
    
    pdb_files = download_all_structures()
    
    if not pdb_files:
        print("✗ Failed to download PDB files. Exiting.")
        return 1
    
    print("="*80)
    print("STEP 3: MEASURE CATALYTIC TRIAD DISTANCES")
    print("="*80 + "\n")
    
    results = analyze_all_structures(pdb_files)
    
    if not results:
        print("✗ Failed to analyze structures. Exiting.")
        return 1
    
    print("="*80)
    print("STEP 4: STRUCTURAL ALIGNMENT (RMSD)")
    print("="*80 + "\n")
    
    rmsd_results = calculate_rmsd_alignments(pdb_files)
    
    print("="*80)
    print("STEP 5: PREPARE FILES FOR CASTp")
    print("="*80 + "\n")
    
    print("Preparing cleaned PDB files for CASTp...")
    for pdb_id, pdb_file in pdb_files.items():
        clean_file = prepare_for_castp(pdb_id, pdb_file)
        print(f"  ✓ {os.path.basename(clean_file)}")
    print()
    
    print("="*80)
    print("STEP 6: GENERATE REPORTS")
    print("="*80 + "\n")
    
    print("Saving analysis results...")
    save_distance_csv(results)
    save_rmsd_csv(rmsd_results)
    save_master_report(results)
    save_markdown_report(results, rmsd_results)
    save_json_report(results, rmsd_results)
    print()
    
    print("="*80)
    print("STEP 7: GENERATE PUBLICATION FIGURES (Optional)")
    print("="*80 + "\n")
    
    generate_pymol_figures()
    
    print("="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80 + "\n")
    
    print("📊 RESULTS SUMMARY:\n")
    print(f"   Directory: {Config.BASE_DIR}/\n")
    print("   Main Files:")
    print("   - catalytic_triad_distances.csv  (CORE DATA)")
    print("   - MASTER_ANALYSIS_REPORT.csv")
    print("   - ANALYSIS_REPORT.md")
    print("   - analysis_results.json")
    if os.path.exists(Config.FIGURES_DIR) and os.listdir(Config.FIGURES_DIR):
        print(f"   - {len(os.listdir(Config.FIGURES_DIR))} publication figures")
    print("\n   Next Steps:")
    print("   1. Review: ANALYSIS_REPORT.md")
    print("   2. CASTp: Upload files from pdb_files/for_castp/")
    print("   3. PDBePISA: Analyze oligomerization interfaces")
    print("\n" + "="*80)
    print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80 + "\n")
    
    return 0

if __name__ == "__main__":
    try:
        exit_code = main()
        sys.exit(exit_code)
    except KeyboardInterrupt:
        print("\n\n⚠ Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n\n✗ ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
