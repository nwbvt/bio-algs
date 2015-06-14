(ns bio-algs.hmm-test
  (:require [clojure.test :refer :all]
            [bio-algs.hmm :refer :all]))

(deftest helper-tests
  (testing "the generation of a transition matrix"
    (is (= (make-matrix
             ["	A	B" "A	0.377	0.623"	"B	0.26	0.74"])
           {\A {\A 0.377 \B 0.623} \B {\A 0.26 \B 0.74}}))))

(deftest hidden-paths
  (testing "finding the probability of a hidden path"
    (is (> 1e-15 (- (p-path "ABABBBAAAA" {\A {\A 0.377 \B 0.623} \B {\A 0.26 \B 0.74}})
                    0.000384928691755))))
  (testing "find the probability of an emission given a path"
    (is (> 1e-15 (- (p-outcome "zzzyxyyzzx" "BAAAAAAAAA" {\A {\x 0.176 \y 0.596 \z 0.228} \B {\x 0.225 \y 0.572 \z 0.203}})
                    3.59748954746e-06))))
  (testing "finding the optimal hidden path using the viterbi algorithm"
    (is (= (optimal-path "xyxzzxyxyy" [\A \B] {\A {\A 0.641 \B 0.359} \B {\A 0.729 \B 0.271}}
                         {\A {\x 0.117 \y 0.691 \z 0.192} \B {\x 0.097 \y 0.42 \z 0.483}})
           "AAABBAAAAA")))
  (testing "find the probability of a given state emitting na particular string"
    (is (> 1e-15 (- (p-hmm-outcome "xzyyzzyzyy" [\A \B] {\A {\A 0.303 \B 0.697} \B {\A 0.831 \B 0.169}}
                                   {\A {\x 0.533 \y 0.065 \z 0.402} \B {\x 0.342 \y 0.334 \z 0.324}})
                    1.1005510319694847e-06)))))

(deftest solving-profile-hmm
  (testing "Solving the profile hmm problem"
    (let [profile (profile-hmm 0.289 ["EBA"
                                      "E-D"
                                      "EB-"
                                      "EED"
                                      "EBD"
                                      "EBE"
                                      "E-D"
                                      "E-D"])]
      (is (= (profile :transition) {:S {:M1 1.0} :M1 {:I1 0.625 :M2 0.375} :I1 {:M2 0.8 :D2 0.2} :M2 {:E 1.0} :D2 {:E 1.0}}))
      (is (= (profile :emission) {:M1 {\E 1.0} :I1 {\B 0.8 \E 0.2} :M2 {\A (double 1/7) \D (double 5/7) \E (double 1/7)}}))
      (is (= (profile :max-state) 2)))))
