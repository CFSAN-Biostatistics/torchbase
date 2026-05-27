version 1.0

task dummy_task {
    input {
        File? contigs
        File? reads
        File? paired1
        File? paired2
        File? interlaced
        File? longreads
    }

    command {
        echo '{"strategy": "balanced", "status": "success"}' > results.json
    }

    output {
        File results = "results.json"
    }
}

workflow balanced_typing {
    input {
        File? contigs
        File? reads
        File? paired1
        File? paired2
        File? interlaced
        File? longreads
    }

    call dummy_task {
        input:
            contigs = contigs,
            reads = reads,
            paired1 = paired1,
            paired2 = paired2,
            interlaced = interlaced,
            longreads = longreads
    }

    output {
        File results = dummy_task.results
    }
}
