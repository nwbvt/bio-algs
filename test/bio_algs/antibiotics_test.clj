(ns bio-algs.antibiotics-test
  (:require [clojure.test :refer :all]
            [bio-algs.antibiotics :refer :all]))

(deftest protein-translation
  (testing "That dna translates into proteins correctly"
    (is (= "MAMAPRTEINSTRING"
           (translate "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA")))))

(deftest peptide-encoding
  (testing "Finding encodings for a peptide"
    (is (= ["ATGGCC" "GGCCAT" "ATGGCC"]
           (find-encoding "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA" "MA")))))

(deftest theoretical-spectrum-test
  (testing "Finding the theoretical spectrum for a peptide"
    (is (= [0 113 114 128 129 227 242 242 257 355 356 370 371 484]
           (theoretical-spectrum "LEQN")))))

(deftest count-peptides-test
  (testing "Counting the number of peptides for a given weight"
    (is (= 3 (count-peptides 128)))))
