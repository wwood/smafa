extern crate assert_cli;
extern crate tempfile;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;
    extern crate tempfile;

    #[test]
    fn test_dna_makedb_and_query() {
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        Assert::main_binary()
            .with_args(&["makedb", "-i", "tests/data/random_3_2.fna", "-d", t])
            .succeeds()
            .unwrap();

        Assert::main_binary()
            .with_args(&["query", "-d", t, "-q", "tests/data/random_3_2.fna"])
            .succeeds()
            .stdout()
            .is("0	0	0	CTT\n\
                1	1	0	AGG\n")
            .unwrap()
    }

    #[test]
    fn test_degenerate_makedb_and_query() {
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        Assert::main_binary()
            .with_args(&["makedb", "-i", "tests/data/degenerate.fna", "-d", t])
            .succeeds()
            .unwrap();

        Assert::main_binary()
            .with_args(&[
                "query",
                "-d",
                t,
                "-q",
                "tests/data/degenerate.fna",
                "--max-divergence",
                "3",
            ])
            .succeeds()
            .stdout()
            .is("0	0	0	CTTNGG\n\
                1	1	0	AGGTGA\n\
                2	2	0	NACTTT\n")
            .unwrap()
    }

    #[test]
    fn test_query_max_divergence_unlimited() {
        Assert::main_binary()
            .with_args(&[
                "query",
                "-d",
                "tests/data/random_3_2.fna.smafadb",
                "-q",
                "tests/data/random_3_2.fna",
                "--max-divergence",
                "99",
                "--max-num-hits",
                "99",
            ])
            .succeeds()
            .stdout()
            .is("0	0	0	CTT\n\
                0	1	3	AGG\n\
                1	1	0	AGG\n\
                1	0	3	CTT\n")
            .unwrap()
    }

    #[test]
    fn test_query_max_divergence_limited() {
        Assert::main_binary()
            .with_args(&[
                "query",
                "-d",
                "tests/data/random_3_2.fna.smafadb",
                "-q",
                "tests/data/random_3_2.fna",
                "--max-divergence",
                "2",
                "--max-num-hits",
                "99",
            ])
            .succeeds()
            .stdout()
            .is("0	0	0	CTT\n\
                1	1	0	AGG\n")
            .unwrap()
    }

    #[test]
    fn test_query_max_divergence_equal() {
        Assert::main_binary()
            .with_args(&[
                "query",
                "-d",
                "tests/data/random_3_2.fna.smafadb",
                "-q",
                "tests/data/random_3_2.fna",
                "--max-divergence",
                "3",
                "--max-num-hits",
                "99",
            ])
            .succeeds()
            .stdout()
            .is("0	0	0	CTT\n\
                0	1	3	AGG\n\
                1	1	0	AGG\n\
                1	0	3	CTT\n")
            .unwrap()
    }

    #[test]
    fn test_query_max_num_hits1() {
        Assert::main_binary()
            .with_args(&[
                "query",
                "-d",
                "tests/data/random_3_2.fna.smafadb",
                "-q",
                "tests/data/random_3_2.fna",
                "--max-num-hits",
                "1",
            ])
            .succeeds()
            .stdout()
            .is("0	0	0	CTT\n\
                1	1	0	AGG\n")
            .unwrap()
    }

    #[test]
    fn test_query_max_num_hits_more() {
        Assert::main_binary()
            .with_args(&[
                "query",
                "-d",
                "tests/data/random_3_2.fna.smafadb",
                "-q",
                "tests/data/random_3_2.fna",
                "--max-num-hits",
                "99",
            ])
            .succeeds()
            .stdout()
            .is("0	0	0	CTT\n\
                0	1	3	AGG\n\
                1	1	0	AGG\n\
                1	0	3	CTT\n")
            .unwrap()
    }

    #[test]
    fn test_fna_count() {
        Assert::main_binary()
            .with_args(&["count", "-i", "tests/data/random_3_2.fna"])
            .succeeds()
            .stdout()
            .is("[{\"path\":\"tests/data/random_3_2.fna\",\"num_reads\":2,\"num_bases\":6}]\n")
            .unwrap()
    }

    #[test]
    fn test_fq_gz_count() {
        Assert::main_binary()
            .with_args(&["count", "-i", "tests/data/random_30_4.fq.gz"])
            .succeeds()
            .stdout()
            .is("[{\"path\":\"tests/data/random_30_4.fq.gz\",\"num_reads\":4,\"num_bases\":120}]\n")
            .unwrap()
    }

    #[test]
    fn test_limit_per_sequence_max_num_hits_2_no_limit() {
        Assert::main_binary()
            .with_args(&[
                "query",
                "-d",
                "tests/data/random_3_2_one_repeated.fna.smafadb",
                "-q",
                "tests/data/random_3_2.fna",
                "--max-num-hits",
                "99",
            ])
            .succeeds()
            .stdout()
            .is("0	0	0	CTT\n\
                0	1	3	AGG\n\
                0	2	3	AGG\n\
                1	1	0	AGG\n\
                1	2	0	AGG\n\
                1	0	3	CTT\n")
            .unwrap()
    }

    #[test]
    fn test_limit_per_sequence_max_num_hits_2_limit1() {
        Assert::main_binary()
            .with_args(&[
                "query",
                "-d",
                "tests/data/random_3_2_one_repeated.fna.smafadb",
                "-q",
                "tests/data/random_3_2.fna",
                "--max-num-hits",
                "99",
                "--limit-per-sequence",
                "1",
            ])
            .succeeds()
            .stdout()
            .is("0	0	0	CTT\n\
                0	1	3	AGG\n\
                1	1	0	AGG\n\
                1	0	3	CTT\n")
            .unwrap()
    }

    // #[test]
    // fn test_db_version_incompatibility(){
    //     Assert::main_binary()
    //         .with_args(&[
    //             "query",
    //             "tests/data/singlem_plot_test.fna_version1_db",
    //             "tests/data/singlem_plot_test.fna"])
    //         .fails()
    //         .stderr().contains(
    //             "Failed to deserialize database, perhaps due to an incompatible database type")
    //         .unwrap();
    // }

    // #[test]
    // fn test_protein_makedb_and_query(){
    //     let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
    //     let t = tf.path().to_str().unwrap();
    //     Assert::main_binary()
    //         .with_args(&[
    //             "makedb",
    //             "tests/data/hello_world_amino_acids.fna",
    //             t,
    //             "--amino-acid"
    //         ]).succeeds().unwrap();

    //     Assert::main_binary()
    //         .with_args(&[
    //             "query",
    //             t,
    //             "tests/data/hello_world_amino_acids.fna",
    //             "-d",
    //             "2"
    //         ]).succeeds()
    //         .stdout().is("1	VVVVVRTWCHHHHH	VVVVVRTKGHHHHH	2	14\n\
    //                       1	VVVVVRTWCHHHHH	VVVVVRTWCHHHHH	0	14\n\
    //                       2	VVVVVRTKGHHHHH	VVVVVRTKGHHHHH	0	14\n\
    //                       2	VVVVVRTKGHHHHH	VVVVVRTWCHHHHH	2	14\n\
    //                       3	VVVVVRAAAHHHHH	VVVVVRAAAHHHHH	0	14\n").unwrap()
    // }

    // #[test]
    // fn test_protein_makedb_and_query_shorter(){
    //     let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
    //     let t = tf.path().to_str().unwrap();
    //     Assert::main_binary()
    //         .with_args(&[
    //             "makedb",
    //             "tests/data/four_amino_acids.faa",
    //             t,
    //             "--amino-acid"
    //         ]).succeeds().unwrap();

    //     Assert::main_binary()
    //         .with_args(&[
    //             "query",
    //             t,
    //             "tests/data/four_amino_acids.faa",
    //             "-d",
    //             "2",
    //             "-k",
    //             "2",
    //         ]).succeeds()
    //         .stdout().is("1	RTWC	RTKG	2	4\n\
    //                       1	RTWC	RTWC	0	4\n\
    //                       2	RTKG	RTKG	0	4\n\
    //                       2	RTKG	RTWC	2	4\n\
    //                       3	RAAA	RAAA	0	4\n").unwrap()
    // }

    // #[test]
    // fn test_translated_makedb_and_query_shorter(){
    //     let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
    //     let t = tf.path().to_str().unwrap();
    //     Assert::main_binary()
    //         .with_args(&[
    //             "makedb",
    //             "tests/data/to_translate.fna",
    //             t,
    //             "--translate"
    //         ]).succeeds().unwrap();

    //     Assert::main_binary()
    //         .with_args(&[
    //             "query",
    //             t,
    //             "tests/data/to_translate.fna",
    //             "--translate",
    //             "-d",
    //             "0",
    //             "-k",
    //             "2",
    //         ]).succeeds()
    //         .stdout().is("1	NCE	NCE	0	3\n\
    //                       1	NCE	NCE	0	3\n\
    //                       2	NCE	NCE	0	3\n\
    //                       2	NCE	NCE	0	3\n").unwrap()
    // }

    // #[test]
    // fn test_cluster_amino_acids() {
    //     let fasta = "tests/data/hello_world_amino_acids.fna";
    //     let mut res = vec!();
    //     smafa::cluster(fasta, 2, &mut res, smafa::SequenceInputType::Translated, 5);
    //     let expected = "C	0	2	*	*	*	*	*	1	*\n\
    //                     C	1	1	*	*	*	*	*	3	*\n\
    //                     H	0	14	85.7	*	*	*	14M	2	1\n";
    //     assert_eq!(expected, String::from_utf8(res).unwrap());
    // }
}
