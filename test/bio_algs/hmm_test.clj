(ns bio-algs.hmm-test
  (:require [clojure.test :refer :all]
            [bio-algs.hmm :refer :all]))

(deftest hidden-paths
  (is (> 1e-15 (- (p-path "ABABBBAAAA" {\A {\A 0.377 \B 0.623} \B {\A 0.26 \B 0.74}})
                  0.000384928691755))))
