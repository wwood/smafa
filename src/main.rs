extern crate bio;
#[macro_use]
extern crate serde_derive;
extern crate bincode;
extern crate serde;
extern crate serde_json;

use std::str;

use bio::alphabets::dna;
use bio::data_structures::fmindex::*;
use bio::data_structures::suffix_array::{suffix_array, SuffixArray};
use bio::data_structures::bwt::*;
use bio::io::fasta;

use std::env;
use std::io;
use std::io::prelude::*;
use std::io::BufWriter;
use std::fs::File;
use std::time::Instant;
use std::rc::Rc;
use std::collections::HashSet;

fn main(){
    let args: Vec<String> = env::args().collect();

    match args.len() {
        2 => {
            match &args[1] as &str {
                "makedb" => makedb(),
                _ => print_help()
            }
        },
        4 => {
            match &args[1] as &str {
                "test2" => test2(&args[2], &args[3]),
                //"query" => query(&args[2]),
                _ => print_help()
            }
        },
        _ => panic!("Unexpected number of arguments")
    }
}

fn print_help(){
    panic!("Not implemented!");
    // println!("Usage:");
    // println!("");
    // println!(" smafa makedb <reference_fasta>");
    // println!("      or");
    // println!(" genome_assigner stage1 <sorted_indexed_bam>");
    // println!("      or");
    // println!(" genome_assigner missing <sorted_indexed_bam> <reference_fasta> <excluded_positions>");
    // println!("");
}

fn makedb(){
    let mut text = vec![];

    let reader = fasta::Reader::new(io::stdin());

    let mut now = Instant::now();
    let mut new_now;

    {
        let mut i: u64 = 0;
        for record in reader.records() {
            text.extend(record.unwrap().seq().iter().cloned());
            text.extend_from_slice(b"$");
            i+=1;
            if i % 10000 == 0 {
                print!(".")
            }
        }
        println!();
    }

    // println!("total: {:?}", text)
    let alphabet = dna::n_alphabet();
    new_now = Instant::now(); println!("reading {:?}", new_now.duration_since(now)); now = new_now;
    let sa = suffix_array(&text);
    new_now = Instant::now(); println!("suffix array {:?}", new_now.duration_since(now)); now = new_now;
    let bwt = bwt(&text, &sa);
    new_now = Instant::now(); println!("bwt {:?}", new_now.duration_since(now)); now = new_now;
    let less = less(&bwt, &alphabet);
    new_now = Instant::now(); println!("less {:?}", new_now.duration_since(now)); now = new_now;
    let occ = Occ::new(&bwt, 3, &alphabet);
    new_now = Instant::now(); println!("occ {:?}", new_now.duration_since(now)); now = new_now;
    let fm = FMIndex::new(&bwt, &less, &occ);
    new_now = Instant::now(); println!("fmindex {:?}", new_now.duration_since(now)); //now = new_now;

    //let serialize_me = SaveableFMIndex;

    {
        let serialized_fm = serde_json::to_string(&fm).unwrap();
        let f = File::create("db.fm").unwrap();
        let mut writer = BufWriter::new(f);
        writer.write(serialized_fm.as_bytes()).unwrap();
    }
    {
        let serialized_sa = serde_json::to_string(&sa).unwrap();
        let f = File::create("db.sa").unwrap();
        let mut writer = BufWriter::new(f);
        writer.write(serialized_sa.as_bytes()).unwrap();
    }
}

// #[derive(Serialize, Deserialize)]
// pub struct Example {
//     fmindex: FMIndex<Rc<BWT>, Rc<Less>, Rc<Occ>>
// }

// fn query(query_sequence: &str){
//     // read in db files
//     let mut fmindex_buf = vec!();
//     //let fmindex: FMIndex<DerefBWT<Target=BWT>, DerefLess<Target=Less>, DerefOcc<Target=Occ>>;
//     //let fmindex: FMIndex<DerefBWT<Target=BWT>, Less, Occ>;
//     //let fmindex: <DerefBWT<Target=BWT>, DerefLess<Target=Less>, DerefOcc<Target=Occ>>;
//     //suffix_array;
//     {
//         let mut f = File::open("db.fm").unwrap();
//         f.read_to_end(&mut fmindex_buf).unwrap();
//         //let inter = ;
//         //println!("inter: {:}", inter);
//     }
//     let fmindex: FMIndex<BWT, Less, Occ> = serde_json::from_str(str::from_utf8(&fmindex_buf).unwrap()).unwrap();
//     //println!("fm: {:?}", fmindex)

//     // query
//     // print
// }

fn test2(db_fasta: &str, query_fasta: &str){
    let mut text = vec![];

    let reader = fasta::Reader::new(File::open(db_fasta).unwrap());

    let mut now = Instant::now();
    let mut new_now;

    {
        let mut i: u64 = 0;
        for record in reader.records() {
            text.extend(record.unwrap().seq().iter().cloned());
            text.extend_from_slice(b"$");
            i+=1;
            if i % 10000 == 0 {
                print!(".")
            }
        }
        println!();
    }

    // println!("total: {:?}", text)
    let alphabet = dna::n_alphabet();
    new_now = Instant::now(); println!("reading {:?}", new_now.duration_since(now)); now = new_now;
    let sa = suffix_array(&text);
    new_now = Instant::now(); println!("suffix array {:?}", new_now.duration_since(now)); now = new_now;
    let bwt = bwt(&text, &sa);
    new_now = Instant::now(); println!("bwt {:?}", new_now.duration_since(now)); now = new_now;
    let less = less(&bwt, &alphabet);
    new_now = Instant::now(); println!("less {:?}", new_now.duration_since(now)); now = new_now;
    let occ = Occ::new(&bwt, 3, &alphabet);
    new_now = Instant::now(); println!("occ {:?}", new_now.duration_since(now)); now = new_now;
    let fm = FMIndex::new(&bwt, &less, &occ);
    new_now = Instant::now(); println!("fmindex {:?}", new_now.duration_since(now)); now = new_now;


    //let serialized = serde_json::to_string(&fm).unwrap();
    //println!("serialized: {:}", serialized);

    let reader = fasta::Reader::new(File::open(query_fasta).unwrap());
    for record in reader.records() {
        let seq = record.unwrap();
        let pattern = seq.seq();
        let mut printed_seqs: HashSet<usize> = HashSet::new();
        for i in 0..11 {
            let intervals = fm.backward_search(pattern[(5*i)..(5*i+5)].iter());
            //new_now = Instant::now(); println!("search {:?}", new_now.duration_since(now)); now = new_now;
            //println!("occ2..");
            // println!("fwd: {:?}", intervals);
            let some = intervals.occ(&sa);
            //new_now = Instant::now(); println!("occ2 {:?}", new_now.duration_since(now));// now = new_now;
            //println!("found: {:} hits", some.len());

            // The text if of length 60, plus the sentinal.
            // So the sequence is located in the text at text[pos/61 .. pos+60] I think.
            for pos in some {
                //println!("pos {:?}", pos);
                let start2 = pos / 61;
                let start = start2 * 61;
                if !printed_seqs.contains(&start) {
                    printed_seqs.insert(start);
                    if pos-start == 0 {
                        //println!("start {:?}", start);
                        let subject = &text[(start..(start+60))];
                        //println!("text {:?}", subject);
                        let mut divergence = 0;
                        let mut total = 0;
                        for i in 0..60 {
                            if subject[i] != pattern[i] {
                                divergence = divergence + 1;
                                if divergence > 5 {
                                    break
                                }
                            }
                            total = total + 1;
                        }
                        if divergence <= 5 {
                            println!("{}\t{}\t{}\t{}\t{}",
                                     seq.id(),
                                     str::from_utf8(pattern).unwrap(),
                                     str::from_utf8(subject).unwrap(),
                                     divergence, total)
                        }
                    }
                }
            }
            //println!()
        }
    }
}

// fn main2(){

//     // let text = b"ATCCC$GGGAT$AGGGA$";
//     // let alphabet = dna::n_alphabet();
//     // let sa = suffix_array(text);
//     // let bwt = bwt(text, &sa);
//     // let less = less(&bwt, &alphabet);
//     // let occ = Occ::new(&bwt, 3, &alphabet);
//     // let fm = FMIndex::new(&bwt, &less, &occ);
//     // //let fmdindex = FMDIndex::from(fm);

//     // let pattern = b"GGGA";
//     // //let intervals = fmdindex.smems(pattern, 2);
//     // let intervals = fm.backward_search(pattern.iter());
//     // println!("fwd: {:?}", intervals);
//     // let positions = intervals.occ(&sa);
//     // println!("fwd: {:?}", positions);
//     // // let forward_positions = intervals[0].forward().occ(&sa);

//     // println!("fwd: {:?}", forward_positions);


//     // let text = b"ATCCC$GGGAT$AGGGA$";
//     // let alphabet = dna::n_alphabet();
//     // let sa = suffix_array(text);
//     // let bwt = bwt(text, &sa);
//     // let less = less(&bwt, &alphabet);
//     // let occ = Occ::new(&bwt, 3, &alphabet);
//     // let fm = FMIndex::new(&bwt, &less, &occ);
//     // search(fm, sa, b"GG")
// }

// pub fn search(fm_index: &FMIndex<DerefBWT<Target=BWT>, DerefLess<Target=Less>, DerefOcc<Target=Occ>>, sa: &SuffixArray, search_pattern: [u8]){
//     let intervals = fm_index.backward_search(search_pattern.iter());
//     println!("fwd: {:?}", intervals.occ(sa));
// }
