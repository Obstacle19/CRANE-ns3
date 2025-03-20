#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 06:46:30 2021

@author: vamsi
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import pylab
from matplotlib.lines import Line2D

results = "./results_myfirst_ss/"
plots_dir = "./plot_myfirst_ss/"
os.makedirs(plots_dir, exist_ok=True)
plt.rcParams.update({'font.size': 18})  

algs = ["dcqcn", "powerInt", "hpcc", "powerDelay", "timely", "dctcp"]
algnames = {
    "dcqcn": "DCQCN",
    "powerInt": "PowerTCP",
    "hpcc": "HPCC",
    "powerDelay": r'$\theta-PowerTCP$',
    "timely": "TIMELY",
    "DCTCP": "DCTCP"
}

#%%

######### BURST ###############

plt.rcParams.update({'font.size': 22})  # 全局字体大小

figlegend = pylab.figure(figsize=(11.5, 1.5))  # 图例
legend_elements = []

# 定义颜色和标签
colorsBurst = ["#1979a9", "red"]
labels = ['Throughput', 'Qlen']

# 创建图例元素
for i in range(2):
    legend_elements.append(Line2D([0], [0], color=colorsBurst[i], lw=6, label=labels[i]))

for alg in algs:
    try:
        # 读取数据
        df = pd.read_csv(results + 'result-' + alg + '.burst', delimiter=' ', usecols=[5, 9, 11, 13],
                         names=["th", "qlen", "time", "power"])

        # 确保数据是一维的
        time = df["time"].values.flatten()
        th = df["th"].values.flatten()
        qlen = df["qlen"].values.flatten()

        # 创建绘图
        fig, ax = plt.subplots(figsize=(15, 6))
        ax.xaxis.grid(True, ls='--')
        ax.yaxis.grid(True, ls='--')

        ax1 = ax.twinx()
        ax.set_yticks([0, 25e9, 50e9, 75e9, 100e9])
        ax.set_yticklabels(["0Gbps", "25Gbps", "50Gbps", "75Gbps", "100Gbps"])
        ax.set_ylabel("Throughput (Gbps)")

        start = 0.13
        xticks = [i * 0.01 + start for i in range(8)]
        ax.set_xticks(xticks)
        ax.set_xticklabels([str(i * 10) for i in range(8)])
        ax.set_xlabel("Time (ms)")
        ax.set_xlim(0.12, 0.20)

        # 绘制 Throughput 和 Queue length
        line1, = ax.plot(time, th, label="Throughput", c='#1979a9', lw=2)  # 吞吐率 
        ax1.set_ylim(0, 1600)
        ax1.set_ylabel("Queue length (KB)")
        line2, = ax1.plot(time, qlen / 1000, c='r', label="Qlen", lw=2)  # 队列长度

        # # 添加图例
        # ax.legend(handles=[line1], loc='upper left')  # 吞吐率的图例
        # ax1.legend(handles=[line2], loc='upper right')  # 队列长度的图例

        # 保存图形
        fig.tight_layout()
        fig.savefig(plots_dir + alg + '.pdf')
        fig.savefig(plots_dir + alg + '.png')

    except Exception as e:
        print(f"Error processing {alg}: {e}")

# 添加图例
figlegend.tight_layout()
figlegend.legend(handles=legend_elements, loc=9, ncol=2, framealpha=0, fontsize=48)