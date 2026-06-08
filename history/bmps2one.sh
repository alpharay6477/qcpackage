#!/bin/bash

# 版权信息
# bmps2one.sh - 将多张图片按指定透明度叠加为一张BMP图片
# 版本: 1.0
# 改进: 添加严格模式、安全临时文件处理、输入验证等

# 设置严格模式
set -euo pipefail

# 检查ImageMagick是否已安装
if ! command -v convert &> /dev/null
then
    echo "错误: ImageMagick未安装，请先安装ImageMagick。"
    exit 1
fi

# 确保至少有一张图片和对应的透明度值
if [ "$#" -lt 2 ]; then
    echo "用法: $0 [选项] 图片文件1 透明度1 [图片文件2 透明度2 ...]"
    echo "选项:"
    echo "  -o, --output <文件名>  指定输出文件名 (默认: output.bmp)"
    exit 1
fi

# 初始化变量
output="output.bmp"  # 默认输出文件名

# 处理命令行选项
while [[ "$1" == -* ]]; do
    case "$1" in
        -o|--output)
            if [ -n "$2" ]; then
                output="$2"
                shift 2
            else
                echo "错误: -o选项需要指定输出文件名"
                exit 1
            fi
            ;;
        *)
            echo "错误: 未知选项 '$1'"
            exit 1
            ;;
    esac
done

# 确保处理选项后至少还有一张图片和对应的透明度值
if [ "$#" -lt 2 ]; then
    echo "错误: 至少需要指定一张图片和对应的透明度值"
    echo "用法: $0 [选项] 图片文件1 透明度1 [图片文件2 透明度2 ...]"
    exit 1
fi

# 创建安全的临时文件
temp_file=$(mktemp "temp_XXXXXX.png")

# 清理函数，确保临时文件被删除
cleanup() {
    if [ -f "$temp_file" ]; then
        rm -f "$temp_file"
    fi
}

# 注册清理函数，确保脚本退出时执行
trap cleanup EXIT

# 处理第一张图片
first_image=$1
first_opacity=$2

# 验证透明度值是否为有效的数字
if ! [[ "$first_opacity" =~ ^[0-9]+$ ]]; then
    echo "错误: 透明度值必须是非负整数，当前值: '$first_opacity'"
    exit 1
fi

# 验证图片文件是否存在
if [ ! -f "$first_image" ]; then
    echo "错误: 图片文件不存在: '$first_image'"
    exit 1
fi

shift 2

# 将第一张图片的透明度设置为指定值
echo "正在处理第一张图片: $first_image (透明度: $first_opacity%)"
convert "$first_image" -alpha set -channel A -evaluate set "$first_opacity"% "$temp_file"

# 循环处理剩余的图片和透明度
counter=2
while [ "$#" -gt 0 ]; do
    image=$1
    opacity=$2
    
    # 验证透明度值是否为有效的数字
    if ! [[ "$opacity" =~ ^[0-9]+$ ]]; then
        echo "错误: 透明度值必须是非负整数，当前值: '$opacity'"
        exit 1
    fi
    
    # 验证图片文件是否存在
    if [ ! -f "$image" ]; then
        echo "错误: 图片文件不存在: '$image'"
        exit 1
    fi
    
    echo "正在处理第 $counter 张图片: $image (透明度: $opacity%)"
    
    # 将当前图片的透明度设置为指定值
    convert "$image" -alpha set -channel A -evaluate set "$opacity"% miff:- |\
    # 叠加当前图片到临时文件上
    convert "$temp_file" miff:- -gravity center -compose over -composite "$temp_file"
    
    shift 2
    counter=$((counter + 1))
done

# 将最终的叠加结果保存为输出文件
echo "正在保存结果到: $output"
mv "$temp_file" "$output"

echo "成功: 输出图片已保存到: $output"


