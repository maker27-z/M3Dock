# M3Dock
# 程序运行与对接准备说明

## 一、运行环境配置

在服务器上安装 **Boost 1.70.0** 库：

- 安装教程参考：
  👉 https://blog.csdn.net/weixin_40998263/article/details/144429066

-  注意：
  - 建议使用管理员权限执行安装：
    ```bash
    sudo ...
    ```
  - 否则可能导致安装失败

---

## 二、对接前的准备

### 蛋白质相关文件
需要准备以下文件：

- `protein.pdb`
- `protein.pdbqt`

---

### 配体相关文件
需要准备以下文件：

- `ligand.mol`
- `ligand.pdbqt`

---

### VDW 能量评分函数所需文件

为了计算 **vdw 能量评分函数**，需要通过 **AutoDock4** 生成相关文件：

- 教程参考：
  👉 https://autodock-vina.readthedocs.io/en/latest/docking_basic.html

---

### 配置文件

- `config.txt`

---

### 参数文件目录

以下目录中的文件是必须的：

MogaDock/LSHADE_Adam_final_nsga2/parameter

---

## 三、示例文件位置

示例相关文件可在以下目录中查看：

./cmake-build-debug

---

## 四、执行命令示例

```bash
./cmake-build-debug/LSHADE_Adam_final \
--config ./example/1qf1/41qf1_config.txt \
--out ./out/1qf1/1qf1_ligand_out.pdbqt \
--log ./out/1qf1/1qf1_ligand_out.log
