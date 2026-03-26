#!/bin/bash

# 分子对接任务运行脚本
# 功能：对每个四位字母数字组合的子目录运行20次分子对接任务

# 定义对接程序路径和参数
DOCKING_PROGRAM="./LSHADE_Adam_final"
CONFIG_PARAM="--config"
OUTPUT_PARAM="--out"
LOG_PARAM="--log"

# 获取当前工作目录
BASE_DIR=$(pwd)

# 获取所有符合要求的子目录（格式为：前一部分1-4位字母数字，下划线，后一部分1-4位字母数字）
SUB_DIRS=()
while IFS= read -r -d $'\0' dir; do
    dir_name=$(basename "$dir")
    if [[ $dir_name =~ ^[[:alnum:]]{1,4}_[[:alnum:]]{1,4}$ ]]; then
        SUB_DIRS+=("$dir")
    fi
done < <(find "$BASE_DIR" -maxdepth 1 -type d -print0)

# 如果没有符合条件的子目录则退出
if [ ${#SUB_DIRS[@]} -eq 0 ]; then
    echo "未找到任何符合条件的子目录！"
    echo "要求：子目录名格式为 '前部分_后部分'，前后部分均由1-4位字母或数字组成（例如：1HPO_AD3, 1GPK_HUP, 1G9V_RQ3）"
    exit 1
fi

echo "找到 ${#SUB_DIRS[@]} 个需要处理的子目录"

# 创建运行日志
RUN_LOG="$BASE_DIR/docking_run_log.txt"
echo "开始分子对接任务，时间: $(date '+%Y-%m-%d %H:%M:%S')" > "$RUN_LOG"

# 计数器
TOTAL_TASKS=$((${#SUB_DIRS[@]} * 20))
COMPLETED_TASKS=0
SKIPPED_RUNS=0

# 主循环
for SUB_DIR in "${SUB_DIRS[@]}"; do
    SUB_DIR_NAME=$(basename "$SUB_DIR")
    LIGAND_ID="$SUB_DIR_NAME"

    # 构建相对路径
    REL_CONFIG="$SUB_DIR_NAME/${SUB_DIR_NAME}_config.txt"
    REL_LIGAND="$SUB_DIR_NAME/${SUB_DIR_NAME}_ligand.pdbqt"

    # 检查所需文件是否存在
    if [ ! -f "$BASE_DIR/$REL_LIGAND" ]; then
        MSG="警告: $SUB_DIR_NAME/ 目录中未找到 ${SUB_DIR_NAME}_ligand.pdbqt 文件，跳过"
        echo "$MSG"
        echo "$MSG" >> "$RUN_LOG"
        continue
    fi
    if [ ! -f "$BASE_DIR/$REL_CONFIG" ]; then
        MSG="警告: $SUB_DIR_NAME/ 目录中未找到 ${SUB_DIR_NAME}_config.txt 文件，跳过"
        echo "$MSG"
        echo "$MSG" >> "$RUN_LOG"
        continue
    fi

    # 检查已完成的运行次数
    COMPLETED=()
    for i in {1..20}; do
        if [ -f "$SUB_DIR/${SUB_DIR_NAME}_${i}_out.pdbqt" ] ; then
            COMPLETED+=("$i")
        fi
    done

    RUNS_NEEDED=()
    for i in {1..20}; do
        if ! [[ " ${COMPLETED[*]} " =~ " $i " ]]; then
            RUNS_NEEDED+=("$i")
        fi
    done

    TOTAL_RUNS_NEEDED=${#RUNS_NEEDED[@]}

    if [ $TOTAL_RUNS_NEEDED -eq 0 ]; then
        MSG="配体 $LIGAND_ID ($SUB_DIR_NAME) 已完成全部20次对接，跳过"
        echo "$MSG"
        echo "$MSG" >> "$RUN_LOG"
        continue
    fi

    MSG="开始处理配体: $LIGAND_ID ($SUB_DIR_NAME) - 已完成 ${#COMPLETED[@]}/20，还需运行 $TOTAL_RUNS_NEEDED 次"
    echo -e "\n***************************************************************"
    echo "$MSG"
    echo "***************************************************************"
    echo -e "\n$MSG" >> "$RUN_LOG"

    # 运行所需的对接
    for RUN in "${RUNS_NEEDED[@]}"; do
        # 构建相对路径的输出文件
        REL_OUT="$SUB_DIR_NAME/${SUB_DIR_NAME}_${RUN}_out.pdbqt"
        REL_LOG="$SUB_DIR_NAME/${SUB_DIR_NAME}_${RUN}_out.log"

        # 执行对接命令（使用相对路径）
        CMD=(
            "$DOCKING_PROGRAM"
            "$CONFIG_PARAM" "$REL_CONFIG"
#            "--ligand" "$REL_LIGAND"
            "$OUTPUT_PARAM" "$REL_OUT"
            "$LOG_PARAM" "$REL_LOG"
        )

        TASK_INFO="配体 $LIGAND_ID 运行 $RUN/20"

        # 运行命令并实时显示输出
        echo -e "\n==================================================================="
        echo "启动对接任务: $TASK_INFO"
        echo "命令: ${CMD[*]}"
        echo "==================================================================="

        START_TIME=$(date +%s)

        echo -e "\n[开始对接输出]"
        "${CMD[@]}" 2>&1 | while IFS= read -r line; do
            # 清理ANSI颜色代码
            CLEAN_LINE=$(echo "$line" | sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,2})?)?[mGK]//g")
            echo "$CLEAN_LINE"
        done

        RETURN_CODE=${PIPESTATUS[0]}

        END_TIME=$(date +%s)
        DURATION=$((END_TIME - START_TIME))
        MINS=$((DURATION / 60))
        SECS=$((DURATION % 60))

        echo -e "\n==================================================================="
        echo "对接任务完成: $TASK_INFO"
        echo "返回代码: $RETURN_CODE"
        echo "运行时间: ${MINS}分 ${SECS}秒"
        echo "===================================================================\n"

        # 记录命令执行结果
        if [ $RETURN_CODE -eq 0 ]; then
            RESULT_MSG="$TASK_INFO - 成功完成对接"
            LOG_MSG="$RESULT_MSG"
        else
            RESULT_MSG="$TASK_INFO - 对接失败, 返回码: $RETURN_CODE"
            LOG_MSG="$RESULT_MSG\n命令: ${CMD[*]}"
        fi

        echo "$RESULT_MSG"
        echo -e "$LOG_MSG" >> "$RUN_LOG"

        # 保存额外信息到单独日志文件
        {
            echo ""
            echo "命令: ${CMD[*]}"
            echo "返回码: $RETURN_CODE"
            echo "任务信息: $TASK_INFO"
        } >> "$BASE_DIR/$REL_LOG"

        COMPLETED_TASKS=$((COMPLETED_TASKS + 1))
        PROGRESS=$(echo "scale=1; $COMPLETED_TASKS * 100 / $TOTAL_TASKS" | bc)

        # 更新进度显示
        echo -e "\n==================================================================="
        echo "当前进度: ${PROGRESS}%"
        echo "总任务数: $TOTAL_TASKS, 已完成: $COMPLETED_TASKS, 跳过的: $SKIPPED_RUNS"
        echo "剩余任务: $((TOTAL_TASKS - COMPLETED_TASKS))"
        echo "===================================================================\n"
    done

    SKIPPED_RUNS=$((SKIPPED_RUNS + 20 - TOTAL_RUNS_NEEDED))
done

# 最终总结
END_TIME=$(date +%s)
TOTAL_DURATION=$((END_TIME - START_TIME))
TOTAL_MINS=$((TOTAL_DURATION / 60))
TOTAL_SECS=$((TOTAL_DURATION % 60))
HOURS=$((TOTAL_MINS / 60))
REMAINING_MINS=$((TOTAL_MINS % 60))

SUMMARY=$(
    echo -e "\n所有对接任务完成!"
    echo "总执行时间: ${HOURS}小时 ${REMAINING_MINS}分 ${TOTAL_SECS}秒"
    echo "总任务数: $TOTAL_TASKS, 成功执行: $COMPLETED_TASKS, 跳过已完成的: $SKIPPED_RUNS"
)

echo "$SUMMARY"
{
    echo ""
    echo "$SUMMARY"
    echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
} >> "$RUN_LOG"

exit 0