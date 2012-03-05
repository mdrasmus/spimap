#!/bin/bash
# make gene example


../../../bin/spimap \
    -a 100.nt.align \
    -s ../../../examples/config/fungi.stree \
    -S ../../../examples/config/fungi.smap \
    -p ../../../examples/train/fungi.params \
    -o 100 \
    -D 0.000564 \
    -L 0.003056 \
    -i 100 \
    --quickiter 1000 \
    -V 1 --log -


cat > 100.yeast.tree <<EOF
(
    (
      (
        (
          (
            (
              (
                YER061C:0.065684,
                spar_6281:0.059258
              )n9:0.024319,
              smik_6662:0.103443
            )n8:0.016745,
            sbay_7039:0.089961
          )n7:0.005255,
          (
            smik_6659:0.092338,
            sbay_7037:0.127706
          )n10:0.014401
        )n6:0.180075,
        CAGL0J02970g:0.290991
      )n5:0.095828,
      scas_g715.48:0.348032
    )n4:0.071532,
    (
      kwal_5828:0.302079,
      (
        ADL072C:0.365623,
        KLLA0C08239g:0.460869
      )n3:0.116995
    )n2:0.054779
)n1;
EOF


~/projects/dlcoal/bin/mpr \
    -s ../../../examples/config/fungi.stree \
    -S ../../../examples/config/fungi.smap \
    -I .tree -O .mpr 100.yeast.tree 


tree-relations -S ../../../examples/config/fungi.smap -s ../../../examples/config/fungi.stree -R .mpr.recon -T .tree 100.yeast.tree


#=============================================================================
# not needed

viewtree -g 100.yeast.tree.svg \
    -l 350 -n \
    -S ../../../examples/config/fungi.smap \
    -s ../../../examples/config/fungi.stree \
    100.yeast.tree 

convert 100.yeast.tree.svg 100.yeast.tree.png

