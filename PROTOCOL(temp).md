采样率：48000 Hz
输入格式：恰好 10000 个 ASCII 字符 0/1，读取时忽略换行和空白
载波频率：10000 Hz
调制方式：BPSK
bit 1：+sin(...)
bit 0：-sin(...)
每 bit 采样数：48
传输速率：约 1000 bit/s

每帧有效数据：100 bit
帧数量：100
CRC：8 bit
生成多项式：110100111
处理顺序：MSB first
CRC 初始值、移位方向和余数附加顺序：必须明确统一

前导：
长度 480 samples
前 240 samples：2 kHz 上升到 10 kHz
后 240 samples：10 kHz 下降到 2 kHz

帧结构：
[480-sample preamble]
[100-bit payload]
[8-bit CRC]
[96-sample silence]

数据载波相位：
每帧数据部分从相位 0 开始
帧内保持连续，不在 callback 边界重置

接收输出：
按帧接收顺序拼接 payload
只把 CRC 正确的 payload 写入结果
最终 OUTPUT.txt 恰好包含 10000 个 0/1