#pragma once
// tqdm-format progress bar for both R and Python builds.
//
// Output format (identical to Python's tqdm default):
//   42%|████████████          | 42/100 [00:05<00:07, 8.33it/s]
//
// R builds:   uses REprintf() so output respects R's console/sink system.
// Python builds (NO_RCPP defined): uses std::fputs(stderr).
// No external library required; same rendering code in both environments.

#ifndef NO_RCPP
#  include <Rcpp.h>   // REprintf
#endif

#include <chrono>
#include <cstdio>
#include <string>

class FastcpdProgress {
 public:
  explicit FastcpdProgress(unsigned int total)
      : total_(total),
        current_(0),
        last_pct_(-1),
        start_(std::chrono::steady_clock::now()) {
    render(0);
  }

  void tick() {
    ++current_;
    int const pct = total_ > 0
        ? static_cast<int>(100.0 * current_ / total_) : 100;
    if (pct > last_pct_ || current_ == total_) {
      render(pct);
      last_pct_ = pct;
    }
  }

 private:
  void render(int pct) const {
    using clock = std::chrono::steady_clock;
    double const elapsed =
        std::chrono::duration<double>(clock::now() - start_).count();
    double const rate = (elapsed > 0.0 && current_ > 0)
        ? current_ / elapsed : 0.0;
    double const eta  = (rate > 0.0 && current_ < total_)
        ? (total_ - current_) / rate : 0.0;

    // Bar: kBarWidth cells, █ (U+2588, UTF-8 \xe2\x96\x88) filled, space empty.
    int const filled = pct * kBarWidth / 100;
    std::string bar;
    bar.reserve(static_cast<size_t>(filled) * 3 +
                static_cast<size_t>(kBarWidth - filled));
    for (int i = 0; i < filled; ++i)         bar += "\xe2\x96\x88";
    for (int i = filled; i < kBarWidth; ++i) bar += ' ';

    char et[32], tt[32], rt[32];
    time_str(et, sizeof(et), elapsed);
    time_str(tt, sizeof(tt), current_ == total_ ? 0.0 : eta);
    rate_str(rt, sizeof(rt), rate);

    char prefix[12], suffix[128];
    std::snprintf(prefix, sizeof(prefix), "\r%3d%%|", pct);
    std::snprintf(suffix, sizeof(suffix), "| %u/%u [%s<%s, %s]",
                  current_, total_, et, tt, rt);

    std::string line = std::string(prefix) + bar + suffix;
    emit(line.c_str());
    if (current_ == total_) emit("\n");
  }

  static void time_str(char* buf, int sz, double sec) {
    if (sec < 0.0) sec = 0.0;
    int const s  = static_cast<int>(sec);
    int const m  = s / 60, sm = s % 60;
    int const h  = m / 60, hm = m % 60;
    if (h > 0) std::snprintf(buf, sz, "%d:%02d:%02d", h, hm, sm);
    else        std::snprintf(buf, sz, "%02d:%02d", m, sm);
  }

  static void rate_str(char* buf, int sz, double rate) {
    if      (rate >= 1.0) std::snprintf(buf, sz, "%.2fit/s", rate);
    else if (rate >  0.0) std::snprintf(buf, sz, "%.2fs/it", 1.0 / rate);
    else                  std::snprintf(buf, sz, "?it/s");
  }

  static void emit(char const* s) {
#ifdef NO_RCPP
    std::fputs(s, stderr);
    std::fflush(stderr);
#else
    REprintf("%s", s);
#endif
  }

  static constexpr int kBarWidth = 20;

  unsigned int const total_;
  unsigned int       current_;
  int                last_pct_;
  std::chrono::steady_clock::time_point start_;
};
