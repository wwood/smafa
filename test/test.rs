extern crate tempfile;
extern crate smafa;

#[cfg(test)]
mod tests {
    #[test]
    fn makedb_and_query100seqs() {
        let mut tf: tempfile::NamedTempFile = tempfile::tempfile().unwrap();
        smafa::makedb(tf.path(), "test/data/4.08.ribosomal_protein_L3_rplC.100random.fna");
        smafa::query(tf.path(), "test/data/4.08.ribosomal_protein_L3_rplC.100random.fna", 5);
    }
}
