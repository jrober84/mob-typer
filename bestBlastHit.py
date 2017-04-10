from Bio.Blast import NCBIXML
import pprint
class bestBlastHit:

    def __init__(self, blast_results):
        self.records = NCBIXML.parse(open(blast_results, "r"))

    def getBestHit(self):
        top_hit_id = ''
        top_hsp = None
        top_hit_score = 0
        query_len = 0
        query_id = ''
        all_hits = dict()

        for record in self.records:
            if not record.alignments:
                continue
            for hit in record.alignments:
                for hsp in hit.hsps:
                    if top_hit_id == '' or top_hit_score < hsp.bits:
                        query_id = record.query
                        query_len = record.query_letters
                        top_hit_id = hit.title
                        top_hsp = hsp
                        top_hit_score = hsp.bits
                    if hit.title not in all_hits:
                        all_hits[hit.title] = dict()
                    all_hits[hit.title][record.query] = {
                        'hit_id': hit.title,
                        'query_len': record.query_letters,
                        'hsp': hsp,
                    }

        return {
            'top_hit_id': top_hit_id,
            'query_id': query_id,
            'query_len': query_len,
            'top_hsp': top_hsp,
            'all_hits': all_hits
        }
