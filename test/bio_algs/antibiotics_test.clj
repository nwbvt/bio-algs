(ns bio-algs.antibiotics-test
  (:require [clojure.test :refer :all]
            [bio-algs.antibiotics :refer :all]))

(deftest protein-translation
  (testing "That dna translates into proteins correctly"
    (is (= "MAMAPRTEINSTRING"
           (translate "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA")))))
