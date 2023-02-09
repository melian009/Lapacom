remvalue(all);
load(linearalgebra);

(%i33) f: r*y*((K-y)/K)-d*x-g*x $

(%i34) g: g*x-a*y-s*y $

(%i35) Jfg: jacobian([f,g],[x,y]) $

(%i36) J00: subst([x=0,y=0],Jfg);
                              [ - d      r     ]
(%o36)                        [                ]
                              [  0   (- s) - a ]
(%i37) J10: subst([x=1,y=0],Jfg);
                        [ - d      r - (1 - a) b     ]
(%o37)                  [                            ]
                        [  0   (- s) + (1 - a) b - a ]
(%i38) J11: subst([x=a,y=1-a],Jfg);
              [                      (a + K - 1) r   (1 - a) r ]
              [ (- d) - (1 - a) a b  ------------- - --------- ]
(%o38)        [                            K             K     ]
              [                                                ]
              [     (1 - a) a b              (- s) - a         ]
(%i39) J: Jfg,x=a,y=1-a;
              [                      (a + K - 1) r   (1 - a) r ]
              [ (- d) - (1 - a) a b  ------------- - --------- ]
(%o39)        [                            K             K     ]
              [                                                ]
              [     (1 - a) a b              (- s) - a         ]
(%i40) determinant(J);
(%o40) ((- d) - (1 - a) a b) ((- s) - a)
                                                     (a + K - 1) r   (1 - a) r
                                      - (1 - a) a b (------------- - ---------)
                                                           K             K


