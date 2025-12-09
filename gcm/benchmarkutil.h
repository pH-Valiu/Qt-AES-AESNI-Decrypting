#ifndef BENCHMARKUTIL_H
#define BENCHMARKUTIL_H

#pragma once
#include <QtCore>
#include <algorithm>
#include <cmath>
#include <x86intrin.h>   // for __rdtsc
#include <gcm/gcm.h>    // for GCM_OUT

class BenchmarkUtil {
public:
    struct Stats {
        quint64 min;
        quint64 max;
        double avg;
        quint64 median;
        double stddev;
    };
    // =============================
    // Random 128-bit generator
    // =============================
    static __m128i rand128() {
        return _mm_set_epi64x(QRandomGenerator::global()->generate64(),
                              QRandomGenerator::global()->generate64());
    }
    static QByteArray randomize(const QByteArray &in) {
        QByteArray out = in;
        for (int i = 0; i < out.size(); i++)
            out[i] = QRandomGenerator::global()->bounded(256);
        return out;
    }
    static QByteArray randQByteArray(int max_size, int min_size = 1) {
        if (max_size <= 0) return QByteArray();

        // Determine random size between 1 and max_size
        int size = QRandomGenerator::global()->bounded(min_size, max_size + 1);

        QByteArray out(size, Qt::Uninitialized);
        for (int i = 0; i < size; i++) {
            out[i] = static_cast<char>(QRandomGenerator::global()->bounded(0, 256));
        }
        return out;
    }

    // =============================
    // Public entry point
    // =============================
    using GfMulFunc = __m128i(*)(const __m128i&, const __m128i&);
    using GfMulFunc512 = void(*)(__m512i* a, __m512i* b, __m512i* res);

    static void run(const QString &name,
                    GfMulFunc512 func,
                    __m512i* fixedA,
                    __m512i* fixedB,
                    __m512i* resPtr,
                    int warmupIters = 20000,
                    int measureIters = 100000,
                    int perCallSamples = 200){

        qInfo() << "====================================";
        qInfo() << "Benchmarking:" << name;
        qInfo() << "====================================";

        // Warm-up
        for (int i = 0; i < warmupIters; i++)
            func(fixedA, fixedB, resPtr);

        // Block measurement
        QVector<quint64> blockTimes;
        blockTimes.reserve(20);
        for (int s = 0; s < 20; s++) {
            quint64 start = __rdtsc();
            for (int i = 0; i < measureIters; i++)
                func(fixedA, fixedB, resPtr);
            quint64 end = __rdtsc();
            blockTimes.append((end - start) / measureIters);
        }
        Stats blockStats = computeStats(blockTimes);

        // Per-call measurement
        qInfo()<<"No random input in per-call measurements";
        QVector<quint64> rawCycles;
        QVector<quint64> singleTimes;
        singleTimes.reserve(perCallSamples);
        for (int i = 0; i < perCallSamples; i++) {
            quint64 start = __rdtsc();
            func(fixedA, fixedB, resPtr);
            quint64 end = __rdtsc();
            singleTimes.append(end - start);
        }

        double cpuGHz = calibratedCPUGHz();
        QVector<quint64> rawNs;
        rawNs.reserve(singleTimes.size());
        for (quint64 c : singleTimes) rawNs.append(quint64(c / cpuGHz));

        qInfo() << "--------- Block measurement ---------";
        printStats(blockStats);
        qInfo() << "--------- Per-call timing ---------";
        printStats(computeStats(singleTimes));
        printHistogramRaw(rawNs, "Per-call histogram (ns)");
    }

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

    // =============================
    // QByteArray-specific benchmarking
    // =============================

    using EncryptFunc = GCM_OUT(*)(const QByteArray &key,
                                       const QByteArray &iv,
                                       const QByteArray &aad,
                                       const QByteArray &plaintext);
    using DecryptFunc = QByteArray(*)(const QByteArray &key,
                                       const QByteArray &iv,
                                       const QByteArray &aad,
                                       const QByteArray &ciphertext,
                                       const QByteArray &tag);

    static void runEncryptBenchmark(const QString &name,
                                    EncryptFunc func,
                                    const QByteArray &key,
                                    const QByteArray &iv,
                                    const QByteArray &aad,
                                    const QByteArray &plaintext,
                                    int warmupIters = 2000,
                                    int measureIters = 10000,
                                    int perCallSamples = 200)
    {
        qInfo() << "Benchmarking encrypt:" << name;

        // Warm-up
        for (int i = 0; i < warmupIters; i++)
            func(key, iv, aad, plaintext);

        // Block measurement
        QVector<quint64> blockTimes;
        blockTimes.reserve(20);
        for (int s = 0; s < 20; s++) {
            quint64 start = __rdtsc();
            for (int i = 0; i < measureIters; i++)
                func(key, iv, aad, plaintext);
            quint64 end = __rdtsc();
            blockTimes.append((end - start) / measureIters);
        }
        Stats blockStats = computeStats(blockTimes);

        // Per-call measurement
        QVector<quint64> rawCycles;
        QVector<quint64> singleTimes;
        singleTimes.reserve(perCallSamples);
        for (int i = 0; i < perCallSamples; i++) {
            QByteArray t1 = randomize(key);
            QByteArray t2 = randQByteArray(256, 96);
            QByteArray t3 = randQByteArray(256, 50);
            QByteArray t4 = randQByteArray(512, 128);
            quint64 start = __rdtsc();
            func(t1, t2, t3, t4);
            quint64 end = __rdtsc();
            singleTimes.append(end - start);
        }

        double cpuGHz = calibratedCPUGHz();
        QVector<quint64> rawNs;
        rawNs.reserve(singleTimes.size());
        for (quint64 c : singleTimes) rawNs.append(quint64(c / cpuGHz));

        qInfo() << "--------- Block measurement ---------";
        printStats(blockStats);
        qInfo() << "--------- Per-call timing ---------";
        printStats(computeStats(singleTimes));
        printHistogramRaw(rawNs, "Per-call histogram (ns)");
    }

    static void runDecryptBenchmark(const QString &name,
                                    DecryptFunc func,
                                    const QByteArray &key,
                                    const QByteArray &iv,
                                    const QByteArray &aad,
                                    const QByteArray &ciphertext,
                                    const QByteArray &tag,
                                    int warmupIters = 2000,
                                    int measureIters = 10000,
                                    int perCallSamples = 200)
    {
        qInfo() << "Benchmarking decrypt:" << name;

        // Warm-up
        for (int i = 0; i < warmupIters; i++)
            func(key, iv, aad, ciphertext, tag);

        // Block measurement
        QVector<quint64> blockTimes;
        blockTimes.reserve(20);
        for (int s = 0; s < 20; s++) {
            quint64 start = __rdtsc();
            for (int i = 0; i < measureIters; i++)
                func(key, iv, aad, ciphertext, tag);
            quint64 end = __rdtsc();
            blockTimes.append((end - start) / measureIters);
        }
        Stats blockStats = computeStats(blockTimes);

        // Per-call measurement
        QVector<quint64> rawCycles;
        QVector<quint64> singleTimes;
        singleTimes.reserve(perCallSamples);
        for (int i = 0; i < perCallSamples; i++) {
            QByteArray t1 = randomize(key);
            QByteArray t2 = randQByteArray(256, 96);
            QByteArray t3 = randQByteArray(256, 50);
            QByteArray t4 = randQByteArray(512, 128);
            QByteArray t5 = randomize(tag);
            quint64 start = __rdtsc();
            func(t1, t2, t3, t4, t5);
            quint64 end = __rdtsc();
            singleTimes.append(end - start);
        }

        double cpuGHz = calibratedCPUGHz();
        QVector<quint64> rawNs;
        rawNs.reserve(singleTimes.size());
        for (quint64 c : singleTimes) rawNs.append(quint64(c / cpuGHz));

        qInfo() << "--------- Block measurement ---------";
        printStats(blockStats);
        qInfo() << "--------- Per-call timing ---------";
        printStats(computeStats(singleTimes));
        printHistogramRaw(rawNs, "Per-call histogram (ns)");
    }


private:

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
                                  int bins = 30)
    {
        if (data.isEmpty()) return;

        // 1. Determine min and max dynamically
        quint64 minv = *std::min_element(data.begin(), data.end());
        quint64 maxv = *std::max_element(data.begin(), data.end());

        if (minv == maxv) {
            // All values are identical, show a single bar
            qInfo() << title;
            qInfo() << "[" << minv << "] " << QString(QChar(0x2588)).repeated(40) << " " << data.size();
            return;
        }

        // 2. Compute step for bins
        quint64 range = maxv - minv;
        quint64 step = std::max<quint64>(1, range / bins);

        // 3. Allocate bins
        QVector<int> hist(bins, 0);

        // 4. Fill bins
        for (quint64 v : data) {
            int idx = std::min<int>((v - minv) / step, bins - 1);
            hist[idx]++;
        }

        // 5. Find max count for bar scaling
        int maxCount = *std::max_element(hist.begin(), hist.end());
        int barWidth = 40;

        // 6. Print histogram
        qInfo() << title;
        qInfo().nospace() << "Range: " << minv << " – " << maxv
                          << "   (bin step = " << step << ")";

        for (int i = 0; i < bins; i++) {
            quint64 lo = minv + i * step;
            quint64 hi = (i == bins - 1) ? maxv : lo + step; // last bin includes max

            int count = hist[i];
            int bars = int(double(count) / maxCount * barWidth);

            QString bar(bars, QChar(0x2588));

            qInfo().nospace()
                << "[" << lo << " – " << hi << "] "
                << bar << " " << count;
        }
    }


};

#endif // BENCHMARKUTIL_H
