#ifndef AVS_API_SEQMODEL_H
#define AVS_API_SEQMODEL_H
#ifdef WSIZ
#  undef WSIZ
#endif

#define M4(a) ((a)[0]>(a)[1]?((a)[0]>(a)[2]?((a)[0]>(a)[3]?(a)[0]:(a)[3]):((a)[2]>(a)[3]?(a)[2]:(a)[3])):((a)[1]>(a)[2]?((a)[1]>(a)[3]?(a)[1]:(a)[3]):((a)[2]>(a)[3]?(a)[2]:(a)[3])))
#include <string.h>
template<typename st_t, int NUMBER = 4>
struct SEQ_MODEL {
    enum {
        STEP = sizeof(st_t) == 1 ? 1 : 8
    };
    enum {
        WSIZ = (1 << 8 * sizeof(st_t)) - 2 * STEP
    };

    SEQ_MODEL();

    SEQ_MODEL(int *start);

    void reset();

    void reset(int *start);

    inline void encodeSymbol(RangeCoder *rc, uint sym);
    inline void addSymbol(uint sym);

    inline void updateSymbol(uint sym);

    inline uint decodeSymbol(RangeCoder *rc);

    inline uint getTopSym(void);

    inline uint getSummFreq(void);

    inline uint getsum(uint sym);

protected:
    void rescaleRare();

    st_t Stats[NUMBER];
};

template<typename st_t, int NUMBER>
SEQ_MODEL<st_t, NUMBER>::SEQ_MODEL() {
    reset();
}

template<typename st_t, int NUMBER>
SEQ_MODEL<st_t, NUMBER>::SEQ_MODEL(int *start) {
    for (int i = 0; i < NUMBER; i++) {
        Stats[i] = start[i];
    }
}

template<typename st_t, int NUMBER>
void SEQ_MODEL<st_t, NUMBER>::reset() {
    for (int i = 0; i < NUMBER; i++)
        Stats[i] = 3 * STEP;
}

template<typename st_t, int NUMBER>
void SEQ_MODEL<st_t, NUMBER>::reset(int *start) {
    for (int i = 0; i < NUMBER; i++) {
        Stats[i] = start[i];
    }
}

template<typename st_t, int NUMBER>
void SEQ_MODEL<st_t, NUMBER>::rescaleRare() {
    for (int i = 0; i < NUMBER; i++) {
        Stats[i] -= (Stats[i] >> 1);
    }
}

template<typename st_t, int NUMBER>
uint SEQ_MODEL<st_t, NUMBER>::getsum(uint sym) {
    int sum = 0;
    for (int i = 0; i < sym; i++) {
        sum += Stats[i];
    }
    return sum;
}

template<typename st_t, int NUMBER>
inline void SEQ_MODEL<st_t, NUMBER>::encodeSymbol(RangeCoder *rc, uint sym) {
    uint SummFreq = getsum(NUMBER);
    if (SummFreq >= WSIZ) {
        rescaleRare();
        SummFreq = getsum(NUMBER);
    }

    rc->Encode(getsum(sym), Stats[sym], SummFreq);

    Stats[sym] += STEP;
    return;
}

template<typename st_t, int NUMBER>
inline void SEQ_MODEL<st_t, NUMBER>::addSymbol(uint sym){
    uint SummFreq = getsum(NUMBER);
    if (SummFreq >= WSIZ) {
        rescaleRare();
        SummFreq = getsum(NUMBER);
    }

    Stats[sym] += STEP;
}

template<typename st_t, int NUMBER>
inline void SEQ_MODEL<st_t, NUMBER>::updateSymbol(uint sym) {
    int SummFreq = (Stats[0] + Stats[1]) + (Stats[2] + Stats[3]);
    if (SummFreq >= WSIZ) {
        rescaleRare();
    }

    Stats[sym] += STEP;
}


template<typename st_t, int NUMBER>
inline uint SEQ_MODEL<st_t, NUMBER>::getTopSym(void) {
    return M4(Stats);
}

template<typename st_t, int NUMBER>
inline uint SEQ_MODEL<st_t, NUMBER>::getSummFreq(void) {
    int SummFreq = (Stats[0] + Stats[1]) + (Stats[2] + Stats[3]);
    return SummFreq;
}

template<typename st_t, int NUMBER>
inline uint SEQ_MODEL<st_t, NUMBER>::decodeSymbol(RangeCoder *rc) {
    uint SummFreq = getsum(NUMBER);
    if (SummFreq >= WSIZ) {
        rescaleRare();
        SummFreq = getsum(NUMBER);
    }

    uint count = rc->GetFreq(SummFreq);
    uint HiCount = 0;

    st_t *p = Stats;
    for (int i = 0; i < NUMBER; i++, p++) {
        if ((HiCount += *p) > count) {
            rc->Decode(getsum(i), *p, SummFreq);
            Stats[i] += STEP;
            return i;
        }
    }
    return 0;
}

#endif