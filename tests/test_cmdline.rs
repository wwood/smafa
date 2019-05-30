extern crate assert_cli;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;

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
}
