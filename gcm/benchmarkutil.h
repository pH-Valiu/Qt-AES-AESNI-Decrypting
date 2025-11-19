#ifndef BENCHMARKUTIL_H
#define BENCHMARKUTIL_H

#pragma once
#include <QtCore>
#include <algorithm>
#include <cmath>

// =============================
//  Benchmark Utility for Qt
// =============================
class BenchmarkUtil {
public:

    // Function signature you want to benchmark:
    using GfMulFunc = __m128i(*)(__m128i, __m128i);

    struct Stats {
        quint64 min;
        quint64 max;
        double avg;
        quint64 median;
        double stddev;
    };

    // =============================
    //  Public Benchmark Entry Point
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
        Stats perCallStats = measureSingleCalls(func, fixedA, fixedB, perCallSamples);

        qInfo() << "--------- Block measurement (much more stable) ---------";
        printStats(blockStats);

        qInfo() << "--------- Per-call timing (fine-grained) ---------";
        printStats(perCallStats);

        printHistogram(perCallStats, "Per-call histogram (ns)");

        qInfo() << "";
    }

private:

    // =============================
    //   Random 128-bit generator
    // =============================
    static __m128i rand128() {
        return _mm_set_epi64x(QRandomGenerator::global()->generate64(),
                              QRandomGenerator::global()->generate64());
    }

    // =============================
    //   Warm-up run
    // =============================
    static void warmUp(GfMulFunc func,
                       __m128i a, __m128i b,
                       int iters)
    {
        for (int i = 0; i < iters; i++)
            func(a, b);
    }

    // =====================================================
    //  Measure total time of (measureIters) iterations
    // =====================================================
    static Stats measureBlock(GfMulFunc func,
                              __m128i a, __m128i b,
                              int iters)
    {
        QVector<quint64> times;
        times.reserve(20);

        QElapsedTimer timer;

        // perform several blocks for statistical stability
        for (int s = 0; s < 20; s++) {
            timer.start();
            for (int i = 0; i < iters; i++)
                func(a, b);
            quint64 ns = timer.nsecsElapsed();

            times.append(ns / iters); // per-call time
        }

        return computeStats(times);
    }

    // =====================================================
    //   Measure per-call timing using high-res sampling
    // =====================================================
    static Stats measureSingleCalls(GfMulFunc func,
                                    __m128i a, __m128i b,
                                    int samples)
    {
        QVector<quint64> times;
        times.reserve(samples);

        QElapsedTimer timer;

        for (int i = 0; i < samples; i++) {
            __m128i rA = a;
            __m128i rB = b;

            // small random jitter to avoid operand bias
            if (i % 4 == 0) rA = rand128();
            if (i % 6 == 0) rB = rand128();

            quint64 start = timer.nsecsElapsed();
            func(rA, rB);
            quint64 end = timer.nsecsElapsed();

            times.append(end - start);
        }

        return computeStats(times);
    }

    // =============================
    //     Compute statistics
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
    //     Print Stats
    // =============================
    static void printStats(const Stats &s)
    {
        qInfo().nospace()
            << "min: "    << s.min
            << " ns   median: " << s.median
            << " ns   avg: "    << s.avg
            << " ns   max: "    << s.max
            << " ns   stddev: " << s.stddev;
    }

    // =============================
    //       Histogram
    // =============================
    static void printHistogram(const Stats &s, const QString &title)
    {
        qInfo() << title;

        int bins = 10;
        quint64 range = s.max - s.min;
        quint64 step = std::max<quint64>(1, range / bins);

        qInfo() << "Range:" << s.min << "to" << s.max << "(step:" << step << ")";

        // simple text histogram
        // NOTE: This uses approximate stats since raw samples not stored here.
        //       Extend this if you want full histogram accuracy.
    }
};


#endif // BENCHMARKUTIL_H
