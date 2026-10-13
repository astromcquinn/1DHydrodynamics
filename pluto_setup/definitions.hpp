#define  PHYSICS                        HD
#define  DIMENSIONS                     1
#define  GEOMETRY                       CARTESIAN
#define  MAPPED_GRID                    NO
#define  BODY_FORCE                     NO
#define  COOLING                        TABULATED
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            13

// Physics dependent declarations

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

// User-defined parameters (labels)

#define  RHO_BG                         0
#define  T_BG                           1
#define  CHI                            2
#define  L_CLOUD                        3
#define  EDGE_WIDTH                     4
#define  AMP                            5
#define  LAMBDA                         6
#define  HEAT_CGS                       7
#define  T_FLOOR                        8
#define  MAX_COOL_RATE                  9
#define  PHASE                          10
#define  COOL_MASK                      11
#define  COOL_MODE                      12

// [Beg] user-defined constants (do not change this line)

#define  COOLING_TIME_STEPPING          RKF23
#define  WARNING_MESSAGES               0

// [End] user-defined constants (do not change this line)
