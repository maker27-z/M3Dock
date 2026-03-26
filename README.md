# M3Dock：A Collaborative Evolutionary–Gradient Optimization Framework for Multi-Objective Molecular Docking

## 1. Environment Setup

Install **Boost 1.70.0** on the server:

- Installation tutorial:
  👉 https://blog.csdn.net/weixin_40998263/article/details/144429066

- ⚠️ Note:
  - It is recommended to run commands with administrator privileges:
    ```bash
    sudo ...
    ```
  - Otherwise, the installation may fail

---

## 2. Pre-docking Preparation

### Protein Files
The following files are required:

- `protein.pdb`
- `protein.pdbqt`

---

### Ligand Files
The following files are required:

- `ligand.mol`
- `ligand.pdbqt`

---

### Files for VDW Energy Scoring Function

To calculate the **van der Waals (VDW) energy scoring function**, you need to generate related files using **AutoDock4**:

- Tutorial:
  👉 https://autodock-vina.readthedocs.io/en/latest/docking_basic.html

---

### Configuration File

- `config.txt`

---

### Parameter Directory

The following directory is required:
./cmake-build-debug/parameter


---

## 3. Example File Location

Example files can be found in:


./cmake-build-debug


---

## 4. Example Command

```bash
./cmake-build-debug/LSHADE_Adam_final \
--config ./example/1qf1/41qf1_config.txt \
--out ./out/1qf1/1qf1_ligand_out.pdbqt \
--log ./out/1qf1/1qf1_ligand_out.log
