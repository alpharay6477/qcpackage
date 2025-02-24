#!/bin/bash

# 检查ImageMagick是否已安装
if ! command -v convert &> /dev/null
then
    echo "ImageMagick未安装，请先安装ImageMagick。"
    exit 1
fi

# 确保至少有一张图片和对应的透明度值
if [ "$#" -lt 2 ]; then
    echo "用法: $0 图片文件1 透明度1 [图片文件2 透明度2 ...]"
    exit 1
fi

# 初始化变量
output="output.bmp"  # 输出文件名
temp_file="temp.png"  # 临时文件，用于叠加

# 处理第一张图片
first_image=$1
first_opacity=$2
shift 2

# 将第一张图片的透明度设置为指定值
convert "$first_image" -alpha set -channel A -evaluate set "$first_opacity"% "$temp_file"

# 循环处理剩余的图片和透明度
while [ "$#" -gt 0 ]; do
    image=$1
    opacity=$2
    shift 2

    # 将当前图片的透明度设置为指定值
    convert "$image" -alpha set -channel A -evaluate set "$opacity"% miff:- |\
    # 叠加当前图片到临时文件上
    convert "$temp_file" miff:- -gravity center -compose over -composite "$temp_file"
done

# 将最终的叠加结果保存为输出文件
mv "$temp_file" "$output"

echo "输出图片已保存到: $output"

