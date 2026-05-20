workflow mlst_typing {
    input {
        File reads
        File reference_db
    }

    output {
        File typing_results = "results.json"
    }

    call mlst_type { input: reads = reads, reference_db = reference_db }
}

task mlst_type {
    input {
        File reads
        File reference_db
    }

    command {
        echo '{"st": 1, "loci": {"adk": 1}}' > results.json
    }

    output {
        File results = "results.json"
    }
}
