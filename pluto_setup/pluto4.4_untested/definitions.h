#define  PHYSICS                        HD
#define  DIMENSIONS                     1
#define  COMPONENTS                     1
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        TABULATED
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        1
#define  USER_DEF_PARAMETERS            12

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

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

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY                   CONST_mp
#define  UNIT_LENGTH                    (1.e3*CONST_pc)
#define  UNIT_VELOCITY                  1.e7

/* [End] user-defined constants (do not change this line) */
