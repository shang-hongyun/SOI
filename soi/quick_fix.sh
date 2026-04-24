#!/bin/bash
set -e

# Fix _apply_rearrangements
python3 -c "
import re
with open('evolution_simulator_ak.py') as f:
    s = f.read()

old = '''    def _apply_rearrangements(self, karyo, centros, node_name, branch_length):
        n_inv = poisson_sample(self.rng, self.inv_rate * branch_length)
        n_rt = poisson_sample(self.rng, self.rt_rate * branch_length)
        n_ncf = poisson_sample(self.rng, self.ncf_rate * branch_length)
        n_eej = poisson_sample(self.rng, self.eej_rate * branch_length)
        n_fis = poisson_sample(self.rng, self.fission_rate * branch_length)
        n_gain = poisson_sample(self.rng, self.gene_gain_rate * branch_length)
        n_tdup = poisson_sample(self.rng, self.tandem_dup_rate * branch_length)
        n_ddup = poisson_sample(self.rng, self.dispersed_dup_rate * branch_length)
        n_utrans = poisson_sample(self.rng, self.unidir_trans_rate * branch_length)

        event_list = (
            [\"inv\"] * n_inv + [\"rt\"] * n_rt + [\"ncf\"] * n_ncf +
            [\"eej\"] * n_eej + [\"fis\"] * n_fis + [\"gain\"] * n_gain +
            [\"tdup\"] * n_tdup + [\"ddup\"] * n_ddup + [\"utrans\"] * n_utrans
        )
        if not event_list:
            return
        self.rng.shuffle(event_list)

        applied = defaultdict(int)
        for e in event_list:
            if e == \"inv\" and self._apply_inversion(karyo, node_name):
                applied[\"inv\"] += 1
            elif e == \"rt\" and self._apply_rt(karyo, node_name):
                applied[\"rt\"] += 1
            elif e == \"ncf\" and self._apply_ncf(karyo, node_name):
                applied[\"ncf\"] += 1
            elif e == \"eej\" and self._apply_eej(karyo, node_name):
                applied[\"eej\"] += 1
            elif e == \"fis\" and self._apply_fission(karyo, node_name):
                applied[\"fis\"] += 1
            elif e == \"gain\" and self._apply_gene_gain(karyo, node_name):
                applied[\"gain\"] += 1
            elif e == \"tdup\" and self._apply_tandem_dup(karyo, node_name):
                applied[\"tdup\"] += 1
            elif e == \"ddup\" and self._apply_dispersed_dup(karyo, node_name):
                applied[\"ddup\"] += 1
            elif e == \"utrans\" and self._apply_unidir_trans(karyo, node_name):
                applied[\"utrans\"] += 1

        self.events.append({\"node\": node_name, \"type\": \"rearrangements\",
                            \"branch_length\": branch_length,
                            \"sampled\": {\"inv\": n_inv, \"rt\": n_rt, \"ncf\": n_ncf,
                                        \"eej\": n_eej, \"fis\": n_fis, \"gain\": n_gain,
                                        \"tdup\": n_tdup, \"ddup\": n_ddup,
                                        \"utrans\": n_utrans},
                            \"applied\": dict(applied)})'''

new = '''    def _apply_rearrangements(self, karyo, centros, node_name, branch_length):
        n_inv = poisson_sample(self.rng, self.inv_rate * branch_length)
        n_rt = poisson_sample(self.rng, self.rt_rate * branch_length)
        n_ncf = poisson_sample(self.rng, self.ncf_rate * branch_length)
        n_eej = poisson_sample(self.rng, self.eej_rate * branch_length)
        n_fis = poisson_sample(self.rng, self.fission_rate * branch_length)
        n_gain = poisson_sample(self.rng, self.gene_gain_rate * branch_length)
        n_loss = poisson_sample(self.rng, self.gene_loss_rate * branch_length)
        n_tdup = poisson_sample(self.rng, self.tandem_dup_rate * branch_length)
        n_ddup = poisson_sample(self.rng, self.dispersed_dup_rate * branch_length)
        n_utrans = poisson_sample(self.rng, self.unidir_trans_rate * branch_length)
        n_frac = poisson_sample(self.rng, self.frac_rate * branch_length)
        n_sdel = poisson_sample(self.rng, self.seg_del_rate * branch_length)
        n_sdup = poisson_sample(self.rng, self.seg_dup_rate * branch_length)
        n_cht = poisson_sample(self.rng, self.chromothripsis_rate * branch_length)

        event_list = (
            [\"inv\"] * n_inv + [\"rt\"] * n_rt + [\"ncf\"] * n_ncf +
            [\"eej\"] * n_eej + [\"fis\"] * n_fis + [\"gain\"] * n_gain +
            [\"loss\"] * n_loss + [\"tdup\"] * n_tdup + [\"ddup\"] * n_ddup +
            [\"utrans\"] * n_utrans + [\"frac\"] * n_frac +
            [\"segdel\"] * n_sdel + [\"segdup\"] * n_sdup + [\"cht\"] * n_cht
        )
        if not event_list:
            return
        self.rng.shuffle(event_list)

        applied = defaultdict(int)
        for e in event_list:
            if e == \"inv\" and self._apply_inversion(karyo, node_name, centros):
                applied[\"inv\"] += 1
            elif e == \"rt\" and self._apply_rt(karyo, node_name, centros):
                applied[\"rt\"] += 1
            elif e == \"ncf\" and self._apply_ncf(karyo, node_name, centros):
                applied[\"ncf\"] += 1
            elif e == \"eej\" and self._apply_eej(karyo, node_name, centros):
                applied[\"eej\"] += 1
            elif e == \"fis\" and self._apply_fission(karyo, node_name, centros):
                applied[\"fis\"] += 1
            elif e == \"gain\" and self._apply_gene_gain(karyo, node_name, centros):
                applied[\"gain\"] += 1
            elif e == \"loss\" and self._apply_gene_loss(karyo, node_name, centros):
                applied[\"loss\"] += 1
            elif e == \"tdup\" and self._apply_tandem_dup(karyo, node_name, centros):
                applied[\"tdup\"] += 1
            elif e == \"ddup\" and self._apply_dispersed_dup(karyo, node_name, centros):
                applied[\"ddup\"] += 1
            elif e == \"utrans\" and self._apply_unidir_trans(karyo, node_name, centros):
                applied[\"utrans\"] += 1
            elif e == \"frac\" and self._apply_fractionation(karyo, centros, node_name):
                applied[\"frac\"] += 1
            elif e == \"segdel\" and self._apply_seg_deletion(karyo, node_name, centros):
                applied[\"segdel\"] += 1
            elif e == \"segdup\" and self._apply_seg_duplication(karyo, node_name, centros):
                applied[\"segdup\"] += 1
            elif e == \"cht\" and self._apply_chromothripsis(karyo, node_name, centros):
                applied[\"cht\"] += 1

        self.events.append({\"node\": node_name, \"type\": \"rearrangements\",
                            \"branch_length\": branch_length,
                            \"sampled\": {\"inv\": n_inv, \"rt\": n_rt, \"ncf\": n_ncf,
                                        \"eej\": n_eej, \"fis\": n_fis, \"gain\": n_gain,
                                        \"loss\": n_loss, \"tdup\": n_tdup, \"ddup\": n_ddup,
                                        \"utrans\": n_utrans, \"frac\": n_frac,
                                        \"segdel\": n_sdel, \"segdup\": n_sdup, \"cht\": n_cht},
                            \"applied\": dict(applied)})'''

s = s.replace(old, new)
with open('evolution_simulator_ak.py', 'w') as f:
    f.write(s)
print('Rearrangements fixed')
"
