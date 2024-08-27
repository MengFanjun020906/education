# 为什么需要时频分析
通过以上的两个例子，我们不难发现傅立叶变换的缺陷。

第一个例子（尽管**这两个信号的时域分布完全相反，但是它们的频谱图是完全一致的**。显然，FFT无法捕捉到信号在时域分布上的不同。）告诉我们，傅里叶变换只能获取一段信号总体上包含哪些频率的成分，但是对各成分出现的时刻并无所知。因此时域相差很大的两个信号，可能频谱图一样。

第二个例子（我们能很明显发现在第二个信号中央的部分出现了一个突变扰动。然而在频域图中，这样的变化并没有很好的被捕捉到。注意到红框中部分，显然**傅里叶变换把突变解释为了一系列低成分的高频信号的叠加，并没有很好的反应突变扰动给信号带来的变化**。）告诉我们，对于信号中的突变，傅里叶变换很难及时捕捉。而在有些场合，这样的突变往往是十分重要的。


总而言之，傅里叶变换非常擅长分析那些频率特征均一稳定的平稳信号。但是对于非平稳信号，傅立叶变换只能告诉我们信号当中有哪些频率成分——而这对我们来讲显然是不够的。我们还想知道各个成分出现的时间。知道信号频率随时间变化的情况，各个时刻的瞬时频率及其幅值——这也就是时频分析

所谓时频分析，就是既要考虑到频率特征，又要考虑到时间序列变化。常用的有两种方法：短时傅里叶变化，以及小波变换。

# 短时傅里叶变换 STFT
时傅里叶变换的思路非常直观：既然对整个序列做FFT会丢失时间信息，那我一段一段地做FFT不就行了嘛！这也正是短时傅里叶变换名称的来源，Short Time Fourier Transorm，这里的 Short Time 就是指对一小段序列做 FFT。

那么怎么一段一段处理呢？直接截取信号的一段来做 FFT 吗？一般我们通过加窗的方法来截取信号的片段。定义一个窗函数 $w(t)$，比如这样。
![[Pasted image 20240825222141.png]]
将窗函数位移到某一中心点 $\tau$，再将窗函数和原始信号相乘就可以得到截取后的信号 $y(t)$。
$$y(t)=x(t)⋅w(t−τ)$$
前面提到的直接截取的方法其实就是对信号加一个矩形窗，不过一般我们很少选用矩形窗，因为矩形窗简单粗暴的截断方法会产生的频谱泄露以及[吉布斯现象](https://www.cnblogs.com/smallqing/p/10221784.html)，不利于频谱分析。

对原始信号 $x (t)$  做 STFT 的步骤如下。

首先将将窗口移动到信号的开端位置，此时窗函数的中心位置在$t =\tau_0$
​
 处，对信号加窗处理
 $$y(t)=x(t)⋅w(t−τ_0)$$
 然后进行傅里叶变换
 $$\mathrm{X}(\omega)=\mathcal{F}(\mathrm{y}(\mathrm{t}))=\int_{\infty}^{+\infty}\mathrm{x}(\mathrm{t})\cdot\mathrm{w}(\mathrm{t}-\tau_{0})\mathrm{e}^{-\mathrm{j}\omega\mathrm{t}}\mathrm{d}\mathrm{t}$$
由此得到第一个分段序列的频谱分布 $X(\omega)$。在现实应用中，由于信号是离散的点序列，所以我们得到的是频谱序列 $X[N]$。

为了便于表示，我们在这里定义函数 $S(\omega, \tau)$，它表示，在窗函数中心为 $\tau$ 时，对原函数进行变换后的频谱结果 $X(\omega)$，即：$$\mathrm{S(\omega,\tau)=\mathcal{F}(x(t)\cdot w(t-\tau))=\int_{\infty}^{+\infty}x(t)\cdot w(t-\tau)e^{-j\omega t}dt}$$

对应到离散场景中，$S[\omega, \tau]$ 就是一个二维矩阵，每一列代表了在不同位置对信号加窗，对得到的分段进行傅里叶变换后的结果序列。
![[Pasted image 20240825225320.png]]
完成了对第一个分段的FFT操作后，移动窗函数到 $\tau_1$​。把窗体移动的距离称为 Hop Size。移动距离一般小于窗口的宽度，从而保证前后两个窗口之间存在一定重叠部分，我们管这个重叠叫 Overlap。
重复以上操作，不断滑动窗口、FFT，最终得到从 $\tau_0 \sim \tau_N$上所有分段的频谱结果：
![[Pasted image 20240825225657.png]]
最终我们得到的 $S$，就是 STFT 变换后的结果。

# Matlab算法实现
STFT 的实现如下，算法返回的三个参数：

`f`: m 维向量，表示傅里叶变换后每个点对应的频率值，单位为 Hz
`t`: n 维向量，表示 n 个窗口中心时间$\tau_1 \sim \tau_n$ ​，单位为秒
`STFT`: 一个二维矩阵 [m, n]，每个列向量代表了在对应 $\tau$ 上 FFT 变换的结果
```matlab
function [STFT, f, t] = mystft(x, win, hop, nfft, fs)
	% 计算短时傅里叶变换
    % Input:
    %   x - 一维信号
    %   win - 窗函数
    %   hop - hop size，移动长度
    %   nfft - FFT points
    %   fs - 采样率
    %
    % Output:
    %   STFT - STFT-矩阵 [T, F]
    %   f - 频率向量
    %   t - 时间向量
    
    % 把 x 变为列向量
    x = x(:);
    xlen = length(x);
    wlen = length(win);

    % 窗口数目 L
    L = 1+fix((xlen-wlen)/hop);
    STFT = zeros(nfft, L);
    
    % STFT
    for l = 0:L-1
        % 加窗
        xw = x(1+l*hop : wlen+l*hop).*win;
        
        % FFT计算
        X = fft(xw, nfft);
        X = fftshift(X);

        STFT(:, 1+l) = X(1:nfft);
    end
    
    % 取每个窗口中点的时间点
    t = (wlen/2:hop:wlen/2+(L-1)*hop)/fs;
    %f = (0:nfft-1)*fs/nfft;
    % 频率 (fftshift之后的)
    f = (-nfft/2:nfft/2-1) * (fs/nfft);
    
end
```

## 使用范例
我们这里使用 Case 1 的范例来看看 STFT 效果如何。

为了方便可视化，这里给出了对 STFT 变换后的可视化函数。
```matlab
function PlotSTFT(T,F,S)
    % Plots STFT
    plotOpts = struct();
    plotOpts.isFsnormalized = false;
    plotOpts.cblbl = getString(message('signal:dspdata:dspdata:MagnitudedB'));
    plotOpts.title = 'Short-time Fourier Transform';
    plotOpts.threshold = max(20*log10(abs(S(:))+eps))-60;
    signalwavelet.internal.convenienceplot.plotTFR(T,F,20*log10(abs(S)+eps),plotOpts);
end
```
注意上面这个函数依赖于 signalwavelet 包，如果没有这个依赖的话，我还提供了一个粗糙版的绘图函数，用法基本类似，不过最后一个参数需要把从窗函数对象传进去：
```matlab
function PlotSTFT_2(T, F, S, win)
    wlen = length(win);
    C = sum(win)/wlen;
    S = abs(S)/wlen/C;
    
    S = 20*log10(S + 1e-6);

    %figure(1)
    surf(T, F, S)
    shading interp;
    axis tight;
    view(0, 90);
    %set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
    xlabel('Time, s');
    ylabel('Frequency, Hz');
    title('Amplitude spectrogram of the signal');
    
    hcol = colorbar;
    %set(hcol, 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel(hcol, 'Magnitude, dB');
end

```
对 Case 1 中的两种情况进行分析，代码如下
```matlab
close all; clear; clc;
fs = 1000;
t = 0:1/fs:1 - 1/fs;
% 窗口大小，推荐取 2 的幂次
wlen = 256;
% hop size 即移动步长，一般要取一个小于 wlen 的数，推荐取 2 的幂次
hop = wlen/4;
% FFT 点数，理论上应该不小于wlen，推荐取 2 的幂次
nfft = 256;

x = [10 * cos(2 * pi * 10 * t), 20 * cos(2 * pi * 20 * t),...
        30 * cos(2 * pi * 30 * t), 40 * cos(2 * pi * 40 * t)];
figure;
subplot(2, 2, 1);
plot(x);
% 随便选的一个窗函数
win = blackman(wlen, 'periodic');

[S, f, t] = mystft(x, win, hop, nfft, fs);
subplot(2, 2, 2);
PlotSTFT(t,f,S);

x = x(end:-1:1);
subplot(2, 2, 3);
plot(x);
win = blackman(wlen, 'periodic');

[S, f, t] = mystft(x, win, hop, nfft, fs);
subplot(2, 2, 4);
PlotSTFT(t,f,S);

```
可以看到，在 FFT 中无法区分的频谱图像在 STFT 中区分就非常明显，可以看出按照不同的时间分段，频谱分布的变化。
![[Pasted image 20240826162535.png]]
为了更好地理解，将中间的图做一次三维旋转：
![[Pasted image 20240826161457.png]]
可以非常清晰地看出频率分布随时间的变换。注意到分界线处存在异常的高频成分（就是 STFT 图像中那三条竖线），这是因为时域信号突变导致的高频成分。
如果你仔细分析上面的内容，你会发现短时傅立叶变换也有不容忽视的缺陷。

最明显的一个问题：窗口的宽度该设多少为好呢？为了阐明这个问题的影响，我们做这么一个实验：调整不同 `wlen` 的值，来看看影响。
![[Pasted image 20240826162611.png]]
注意 `S` 的尺寸随 `wlen` 的变换，不难发现一个事实：

**窗太窄，窗内的信号太短，会导致频率分析不够精准，频率分辨率差**，具体表现是黄色的横线越来越宽、越来越模糊
**窗太宽，时域上又不够精细，时间分辨率低**，具体表现是淡蓝色的竖线越来越宽、越来越模糊（还记得吗，竖线表示交界处的突变造成的高频干扰成分）
从定量的角度来看，STFT的时间分辨率取决于滑移宽度 $H$，而频率分辨率则取决于 $\frac{F_s}{H}$ ​ 。显然，一方的增加必然意味着另一方的减小。这就是所谓的时频测不准原理（跟海森堡测不准是一个性质），具体关系为：
$$\Delta\mathrm{t}\cdot\Delta\mathrm{f}\geqslant\frac{1}{4\pi}$$
![[Pasted image 20240826162951.png]]
另外，固定的窗口大小过于死板。对低频信号而言，有可能连一个周期都不能覆盖；对高频信号而言，可能覆盖过多周期，不能反映信号变化。