use std::{io::BufRead, fs::File, error::Error};
use serde::{Serialize, Deserialize};
use std::time::{Instant};
use std::io::{Write, Read};

#[derive(Serialize, Deserialize, Debug)]
struct WindowSet {
    version: u32,
    windows: Vec<Vec<bool>>
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() == 1 {
        return encode();
    } else if args.len() == 2 {
        if args[1] == "decode" {
            return decode(&"windows.cbor");
        } else if args[1] == "query" {
            return query(&"windows.cbor");
        } else if args[1] == "encode" {
            return encode();
        } else {
            panic!("Unknown argument: {}", args[1]);
        }
    } else {
        return Err("Usage: ./main subcommand".into());
    }
}

fn encode() -> Result<(), Box<dyn Error>> {
    // Iterate over lines, creating a vector of bools by one hot encoding
    // the input.
    let encoded = std::io::stdin()
        .lock()
        .lines()
        .map(|line| {
            line.unwrap()
                .chars()
                .map(|c| 
                    match c {
                        'A' => [true,false,false,false,false],
                        'C' => [false,true,false,false,false],
                        'G' => [false,false,true,false,false],
                        'T' => [false,false,false,true,false],
                        _ => [false,false,false,false,true]
                    }
                )
                .flatten()
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let windows = WindowSet { 
        version: 1,
        windows: encoded };

    // Encode
    let mut ferris_file = File::create("windows.cbor")?;
    // open binary file
    ferris_file.write(&postcard::to_allocvec(&windows).unwrap())?;
    // serde_cbor::to_writer(ferris_file, &windows)?;
    return Ok(());
}

fn decode(filename: &str) -> Result<(), Box<dyn Error>> {
    panic!();
    // // Decode
    // let ferris_file = File::open(filename)?;
    // let windows: WindowSet = postcard::from_bytes(ferris_file.read(buf).unwrap())?;

    // // Iterate over windows, decoding them into a string.
    // for window in windows.windows {
    //     let mut s = String::new();
    //     for i in 0..window.len() / 5 {
    //         let slice = &window[i*5..(i+1)*5];
    //         s.push(match slice {
    //             [true, false, false, false, false] => 'A',
    //             [false, true, false, false, false] => 'C',
    //             [false, false, true, false, false] => 'G',
    //             [false, false, false, true, false] => 'T',
    //             _ => 'x'
    //         });
    //     }
    //     println!("{}", s);
    // }
    // return Ok(());
}

fn query(filename: &str) -> Result<(), Box<dyn Error>> {
    // Decode
    let mut start = Instant::now();
    let mut ferris_file = File::open(filename)?;
    let mut buffer = Vec::new();
    ferris_file.read_to_end(&mut buffer)?;
    let windows: WindowSet = postcard::from_bytes(&buffer)?;
    eprintln!("Decoded in {}ms", start.elapsed().as_millis()); start = Instant::now();

    // encode a line from stdin as a vector of bools
    let query_vec = std::io::stdin()
        .lock()
        .lines()
        .next()
        .unwrap()
        .unwrap()
        .chars()
        .map(|c| 
            match c {
                'A' => [true,false,false,false,false],
                'C' => [false,true,false,false,false],
                'G' => [false,false,true,false,false],
                'T' => [false,false,false,true,false],
                _ => [false,false,false,false,true]
            }
        )
        .flatten()
        .collect::<Vec<_>>();

    // Get the minimum distance between the query and each window using xor.
    let distances = windows.windows
        .iter()
        .map(|window| {
            window.iter()
                .zip(query_vec.iter())
                .map(|(a, b)| a ^ b)
                .filter(|x| *x)
                .count()
        })
        .collect::<Vec<_>>();
    eprintln!("Distanced in {}ms", start.elapsed().as_millis()); start = Instant::now();

    // Find the minimum distance.
    let min_distance = distances.iter().min().unwrap();
    eprintln!("Found minimum in {}ms", start.elapsed().as_millis()); start = Instant::now();

    // // Print the windows with the minimum distance.
    // for (i, distance) in distances.iter().enumerate() {
    //     if distance == min_distance {
    //         let mut s = String::new();
    //         for j in 0..windows.windows[i].len() / 5 {
    //             let slice = &windows.windows[i][j*5..(j+1)*5];
    //             s.push(match slice {
    //                 [true, false, false, false, false] => 'A',
    //                 [false, true, false, false, false] => 'C',
    //                 [false, false, true, false, false] => 'G',
    //                 [false, false, false, true, false] => 'T',
    //                 _ => 'x'
    //             });
    //         }
    //         println!("{} {} {}", min_distance, i, s);
    //     }
    // }
    // eprintln!("Printed in {}ms", start.elapsed().as_millis());
    return Ok(());
}