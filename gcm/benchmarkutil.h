#ifndef BENCHMARKUTIL_H
#define BENCHMARKUTIL_H

#pragma once
#include <QtCore>
#include <algorithm>
#include <cmath>
#include <x86intrin.h>   // for __rdtsc

class BenchmarkUtil {
public:

    using GfMulFunc = __m128i(*)(__m128i, __m128i);

    struct Stats {
        quint64 min;
        quint64 max;
        double avg;
        quint64 median;
        double stddev;
    };

    // =============================
    // Public entry point
    // =============================
    static void run(const QString &name,
                    GfMulFunc func,
                    __m128i fixedA,
                    __m128i fixedB,
                    int warmupIters      = 20000,
                    int measureIters     = 100000,
                    int perCallSamples   = 200)
    {
        qInfo() << "====================================";
        qInfo() << "Benchmarking:" << name;
        qInfo() << "====================================";

        warmUp(func, fixedA, fixedB, warmupIters);

        Stats blockStats   = measureBlock(func, fixedA, fixedB, measureIters);

        QVector<quint64> rawCycles;
        Stats perCallStats = measureSingleCallsRDTSC(
            func, fixedA, fixedB,
            perCallSamples, &rawCycles);

        // convert cycles -> ns using approximate CPU freq
        double cpuGHz = calibratedCPUGHz();
        QVector<quint64> rawNs;
        rawNs.reserve(rawCycles.size());
        for (quint64 c : rawCycles)
            rawNs.append(quint64(c / cpuGHz));

        qInfo() << "--------- Block measurement (stable) ---------";
        printStats(blockStats);

        qInfo() << "--------- Per-call timing (RDTSC) ---------";
        printStats(perCallStats);

        printHistogramRaw(rawNs, "Per-call timing histogram (ns)");
    }

private:

    // =============================
    // Random 128-bit generator
    // =============================
    static __m128i rand128() {
        return _mm_set_epi64x(QRandomGenerator::global()->generate64(),
                              QRandomGenerator::global()->generate64());
    }

    // =============================
    // Warm-up run
    // =============================
    static void warmUp(GfMulFunc func,
                       __m128i a, __m128i b,
                       int iters)
    {
        for (int i = 0; i < iters; i++)
            func(a, b);
    }

    // =============================
    // Block measurement using RDTSC
    // =============================
    static Stats measureBlock(GfMulFunc func,
                              __m128i a, __m128i b,
                              int iters)
    {
        QVector<quint64> times;
        times.reserve(20);

        for (int s = 0; s < 20; s++) {
            uint64_t start = __rdtsc();
            for (int i = 0; i < iters; i++)
                func(a, b);
            uint64_t end = __rdtsc();
            times.append((end - start) / iters);
        }

        return computeStats(times);
    }

    // =============================
    // Per-call measurement using RDTSC
    // =============================
    static Stats measureSingleCallsRDTSC(GfMulFunc func,
                                         __m128i a, __m128i b,
                                         int samples,
                                         QVector<quint64>* outRaw)
    {
        QVector<quint64> times;
        times.reserve(samples);

        for (int i = 0; i < samples; i++) {
            __m128i rA = a;
            __m128i rB = b;
            if (i % 4 == 0) rA = rand128();
            if (i % 6 == 0) rB = rand128();

            uint64_t start = __rdtsc();
            func(rA, rB);
            uint64_t end = __rdtsc();

            times.append(end - start);
        }

        if (outRaw) *outRaw = times;
        return computeStats(times);
    }

    // =============================
    // Estimate CPU frequency in GHz
    // =============================
    static double calibratedCPUGHz()
    {
        // Measure wall clock over 100ms
        const int iters = 1000000;
        QElapsedTimer t;
        t.start();
        uint64_t start = __rdtsc();
        for (int i = 0; i < iters; i++) { __asm__ volatile(""); }
        uint64_t end = __rdtsc();
        qint64 ms = t.elapsed();
        if (ms == 0) ms = 1;
        double cyclesPerMs = double(end - start) / ms;
        return cyclesPerMs / 1e6; // GHz
    }

    // =============================
    // Compute statistics
    // =============================
    static Stats computeStats(const QVector<quint64> &v)
    {
        QVector<quint64> x = v;
        std::sort(x.begin(), x.end());

        Stats st;
        st.min = x.first();
        st.max = x.last();

        double sum = 0;
        for (quint64 t : x) sum += t;
        st.avg = sum / x.size();

        st.median = x[x.size() / 2];

        double var = 0;
        for (quint64 t : x) var += (t - st.avg) * (t - st.avg);
        var /= x.size();
        st.stddev = std::sqrt(var);

        return st;
    }

    // =============================
    // Print statistics
    // =============================
    static void printStats(const Stats &s)
    {
        qInfo().nospace()
            << "min: "    << s.min
            << " cycles   median: " << s.median
            << " cycles   avg: "    << s.avg
            << " cycles   max: "    << s.max
            << " cycles   stddev: " << s.stddev;
    }

    // =============================
    // Histogram
    // =============================
    static void printHistogramRaw(const QVector<quint64> &data,
                                  const QString &title,
                                  int bins = 30,
                                  quint64 fixedMin = 5) // <--- new parameter
    {
        if (data.isEmpty()) return;

        quint64 minv = fixedMin;                   // always start at 0ns
        quint64 maxv = *std::max_element(data.begin(), data.end());
        quint64 range = maxv - minv;
        quint64 step  = 1;//std::max<quint64>(1, range / bins);

        QVector<int> hist(bins, 0);
        for (quint64 v : data) {
            int idx = std::min<int>((v - minv) / step, bins - 1);
            hist[idx]++;
        }

        int maxCount = *std::max_element(hist.begin(), hist.end());
        int barWidth = 40;

        qInfo() << title;
        qInfo().nospace() << "Range: " << minv << " – " << maxv
                          << " ns   (bin step = " << step << ")";

        for (int i = 0; i < bins; i++) {
            quint64 lo = minv + i * step;
            quint64 hi = lo + step;

            int count = hist[i];
            int bars = (int)((double)count / maxCount * barWidth);

            QString bar(bars, QChar(0x2588));

            qInfo().nospace()
                << "[" << lo << " – " << hi << " ns] "
                << bar << " " << count;
        }
    }

};

#endif // BENCHMARKUTIL_H
