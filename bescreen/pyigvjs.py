def pyigvjs_one(refgenome,
                refgenemefai,
                bamfile,
                bamfilefai,
                gtffile,
                gtffiletbi,
                locus):
    igvjs_html = f"""\
<div id="igv_div">

    <script type="module">
    
        import igv from "https://cdn.jsdelivr.net/npm/igv@3.0.0/dist/igv.esm.min.js"
    
        const div = document.getElementById("igv_div")
    
        const options =
        {{
        reference:
            {{
            "id": "hg38",
            "fastaURL": "{refgenome}",
            "indexURL": "{refgenemefai}",
            "locus": "{locus}",
            "tracks":
                [
                    {{
                    "type": "annotation",
                    "format": "gtf",
                    "name": "annotations",
                    "url": "{gtffile}",
                    "indexURL": "{gtffiletbi}",
                    "colorBy": "transcript_biotype",
                    "colorTable":
                        {{
                        "protein_coding": "magenta",
                        "*": "lightgrey",
                        }}
                    }},
                    {{
                    "type": "alignment",
                    "format": "bam",
                    "name": "found guides",
                    "showCoverage": false,
                    "url": "{bamfile}",
                    "indexURL": "{bamfilefai}",
                    }}
                ]
            }}
        }}
    
        const browser = await igv.createBrowser(div, options)
    
    </script>
</div>
    """

    return igvjs_html


def pyigvjs_two(refgenome,
                refgenemefai,
                bamfile,
                bamfilefai,
                bamfile_ne,
                bamfilefai_ne,
                gtffile,
                gtffiletbi,
                locus):
    igvjs_html = f"""\
<div id="igv_div">

    <script type="module">
    
        import igv from "https://cdn.jsdelivr.net/npm/igv@3.0.0/dist/igv.esm.min.js"
    
        const div = document.getElementById("igv_div")
    
        const options =
        {{
        reference:
            {{
            "id": "hg38",
            "fastaURL": "{refgenome}",
            "indexURL": "{refgenemefai}",
            "locus": "{locus}",
            "tracks":
                [
                    {{
                    "type": "annotation",
                    "format": "gtf",
                    "name": "annotations",
                    "url": "{gtffile}",
                    "indexURL": "{gtffiletbi}",
                    "colorBy": "transcript_biotype",
                    "colorTable":
                        {{
                        "protein_coding": "magenta",
                        "*": "lightgrey",
                        }}
                    }},
                    {{
                    "type": "alignment",
                    "format": "bam",
                    "name": "editing guides",
                    "showCoverage": false,
                    "url": "{bamfile}",
                    "indexURL": "{bamfilefai}",
                    }},
                    {{
                    "type": "alignment",
                    "format": "bam",
                    "name": "non-editing guides",
                    "showCoverage": false,
                    "url": "{bamfile_ne}",
                    "indexURL": "{bamfilefai_ne}",
                    }}
                ]
            }}
        }}
    
        const browser = await igv.createBrowser(div, options)
    
    </script>
</div>
    """

    return igvjs_html


if __name__ == '__main__':
    pass
