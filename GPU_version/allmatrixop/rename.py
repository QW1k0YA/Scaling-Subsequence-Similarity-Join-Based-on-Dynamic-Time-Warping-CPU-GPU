import os
import logging
from datetime import datetime

# 配置日志系统[1,6](@ref)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('rename.log'),
        logging.StreamHandler()
    ]
)

def rename_cpp_to_cu():
    """将当前目录下所有.cpp文件重命名为.cu文件"""
    current_dir = os.getcwd()
    processed = 0
    skipped = 0
    errors = 0

    for filename in os.listdir(current_dir):
        # 检查.cpp扩展名并过滤目录[3,7](@ref)
        if not filename.endswith('.cpp') or not os.path.isfile(filename):
            continue

        # 分离文件名和扩展名[3,5](@ref)
        base_name = os.path.splitext(filename)[0]
        new_name = f"{base_name}.cu"

        # 构建完整路径[2,6](@ref)
        old_path = os.path.join(current_dir, filename)
        new_path = os.path.join(current_dir, new_name)

        # 检查文件是否已存在[4,7](@ref)
        if os.path.exists(new_path):
            logging.warning(f"跳过已存在文件: {new_name}")
            skipped += 1
            continue

        try:
            # 执行重命名操作[4,8](@ref)
            os.rename(old_path, new_path)
            logging.info(f"重命名成功: {filename} → {new_name}")
            processed += 1
        except Exception as e:
            logging.error(f"重命名失败: {filename} | 错误: {str(e)}")
            errors += 1

    # 生成统计报告[1,6](@ref)
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    report = f"\n=== 执行报告 ({timestamp}) ===\n" \
             f"处理文件: {processed}\n" \
             f"跳过文件: {skipped}\n" \
             f"错误次数: {errors}\n" \
             f"============================"
    logging.info(report)

if __name__ == "__main__":
    rename_cpp_to_cu()