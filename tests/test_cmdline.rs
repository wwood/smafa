extern crate assert_cli;
extern crate tempfile;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;
    extern crate tempfile;

    #[test]
    fn test_db_version_incompatibility(){
        Assert::main_binary()
            .with_args(&[
                "query",
                "tests/data/singlem_plot_test.fna_version1_db",
                "tests/data/singlem_plot_test.fna"])
            .fails()
            .stderr().contains(
                "Failed to deserialize database, perhaps due to an incompatible database type")
            .unwrap();
    }

    #[test]
    fn test_protein_makedb_and_query(){
        let tf: tempfile::NamedTempFile = tempfile::NamedTempFile::new().unwrap();
        let t = tf.path().to_str().unwrap();
        Assert::main_binary()
            .with_args(&[
                "makedb",
                "tests/data/hello_world_amino_acids.fna",
                t,
                "--amino-acid"
            ]).succeeds().unwrap();

        Assert::main_binary()
            .with_args(&[
                "query",
                t,
                "tests/data/hello_world_amino_acids.fna",
                "-d",
                "2"
            ]).succeeds()
            .stdout().is("1	VVVVVRTWCHHHHH	VVVVVRTKGHHHHH	2	14\n\
                          1	VVVVVRTWCHHHHH	VVVVVRTWCHHHHH	0	14\n\
                          2	VVVVVRTKGHHHHH	VVVVVRTKGHHHHH	0	14\n\
                          2	VVVVVRTKGHHHHH	VVVVVRTWCHHHHH	2	14\n\
                          3	VVVVVRAAAHHHHH	VVVVVRAAAHHHHH	0	14\n").unwrap()
    }
}
