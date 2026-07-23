import os
import sys
import subprocess
import time

def print_step(message):
    print(f"\n[>>>] {message}")

def main():
    print("="*60)
    print("   🚀 BIOINFORMATICS CLINICAL SUITE - AUTOMATED PIPELINE")
    print("="*60)

    # 1. Python Sürüm Kontrolü
    print_step("Checking system requirements...")
    if sys.version_info < (3, 8):
        print(" [!] Error: Python 3.8 or higher is required.")
        sys.exit(1)
    else:
        print(f"  - Python Version: {sys.version.split()[0]} [OK]")

    # 2. Gerekli Klasörlerin Kontrolü (Yoksa otomatik oluşturur)
    print_step("Verifying project structure...")
    directories = ["scripts", "results", "results/geo_data"]
    for d in directories:
        os.makedirs(d, exist_ok=True)
        print(f"  - Directory '{d}/' [OK]")

    # 3. MUSCLE Motoru Kontrolü
    muscle_path = os.path.join("scripts", "muscle.exe" if os.name == "nt" else "muscle")
    if not os.path.exists(muscle_path):
        print(f"  [!] Warning: '{muscle_path}' not found. Phylogenetic tree module will use fallback simulations.")
    else:
        print(f"  - MUSCLE Alignment Engine [OK]")

    # 4. Arayüz Dosyasının Kontrolü
    print_step("Preparing to launch the Clinical Decision Support System...")
    interface_script = os.path.join("scripts", "task06_web_interface.py")
    
    if not os.path.exists(interface_script):
        print(f" [!] Error: Interface script '{interface_script}' not found! Please check the path.")
        sys.exit(1)

    print("  - Interface script found [OK]")
    print("\nStarting Streamlit server in 3 seconds... (Press Ctrl+C to stop the system)")
    time.sleep(3)
    
    # 5. Sistemi Ayağa Kaldırma
    try:
        subprocess.run([sys.executable, "-m", "streamlit", "run", interface_script])
    except KeyboardInterrupt:
        print("\n\n[>>>] Pipeline securely terminated by user. Have a great day!")

if __name__ == "__main__":
    main()