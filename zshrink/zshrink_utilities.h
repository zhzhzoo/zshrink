#ifdef DEBUG
#   define WHERESTR  "[%s, l%d]: "
#   define WHEREARG  __FILENAME__, __LINE__
#   define DEBUGPRINT2(...)       do {\
                                    fprintf(stderr, __VA_ARGS__);\
                                  } while(0);
#   define DEBUGPRINT(_fmt, ...)  DEBUGPRINT2(WHERESTR _fmt, WHEREARG, ##__VA_ARGS__)
#else
#   define DEBUGPRINT(...)
#endif
