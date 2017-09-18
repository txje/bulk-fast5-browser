import json
import numpy as np
import h5py
import argparse
import asyncio
import websockets
import traceback

"""
 root
   IntermediateData
     Channel_307
       Meta
         description : Grouper
         range : 2131.14
         digitisation : 8192.0
         offset : 13.0
         sample_rate : 6024.0
         threshold : 0.0
         window : 0
         elimit : 0.0
         smallest_event : 0.0
         scaling_used : 1
       Reads
         Data: {'names':['read_id','read_number','event_index_start','event_index_end','read_start','read_length','classification','pen_classification','modal_classification','current_well_id','local_median','local_median_before','local_median_sd','local_median_dwell','local_range','median','median_before','median_sd','median_dwell','range','drift','flags','read_mean_min','read_mean_max'], 'formats':['S37','<u4','<u8','<u8','<u8','<u4','i1','i1','i1','<u4','<f8','<f8','<f8','<f8','<f8','<f8','<f8','<f8','<f8','<f8','<f8','u1','<f4','<f4'], 'offsets':[0,40,168,176,56,64,280,288,284,292,192,200,208,216,224,232,240,248,256,264,272,296,300,304], 'itemsize':308} (7028,)
       BasecallSummary
         Data: {'names':['read_id','read_number','start_time','duration','basecall_start_index_template','basecall_end_index_template','basecall_start_index_complement','basecall_end_index_complement','basecall_event_start_index_template','basecall_event_end_index_template','basecall_event_start_index_complement','basecall_event_end_index_complement','has_template','has_complement','first_sample_template','duration_template','first_sample_complement','duration_complement','result_segmentation','sequence_length_template','num_events_template','called_events_template','mean_qscore_template','strand_score_template','stay_prob_template','step_prob_template','skip_prob_template','sequence_length_complement','num_events_complement','called_events_complement','mean_qscore_complement','strand_score_complement','stay_prob_complement','step_prob_complement','skip_prob_complement','result_basecall_1d'], 'formats':['S37','<u4','<u8','<u8','<u8','<u8','<u8','<u8','<u8','<u8','<u8','<u8','u1','u1','<u8','<u8','<u8','<u8','i1','<u8','<u8','<u8','<f4','<f4','<f4','<f4','<f4','<u8','<u8','<u8','<f4','<f4','<f4','<f4','<f4','i1'], 'offsets':[4,0,48,56,400,408,432,440,416,424,448,456,101,102,64,72,80,88,100,104,112,120,128,132,136,140,144,152,160,168,176,180,184,188,192,200], 'itemsize':464} (0,)
       Events
         Data: [('start', '<u8'), ('length', '<u4'), ('flags', '<u4'), ('mean', '<f4'), ('variance', '<f4')] (875392,)
       BasecallSequence
         Data: [('base', 'i1'), ('qscore', 'i1')] (0,)
     ...
   StateData
     Channel_307
       States
         Data: [('acquisition_raw_index', '<u8'), ('analysis_raw_index', '<u8'), ('summary_state', '<i4')] (1266,)
       Meta
     ...
   MultiplexData
     Channel_307
       Multiplex
         Data: [('approx_raw_index', '<u8'), ('well_id', '<i4')] (19280,)
       Meta
     ...
   Device
     MetaData
       Data: [('asic_temperature', '<f4'), ('heatsink_temperature', '<f4'), ('bias_voltage', '<i2')] (34894336,)
     AsicCommands
       Data: {'names':['frame_number','command','channel_input_reason'], 'formats':['<u8','V512','V512'], 'offsets':[0,16,528], 'itemsize':1040} (19280,)
   UniqueGlobalKey
     context_tags
       experiment_duration_set : 2880
       experiment_type : genomic_dna
       fast5_output_fastq_in_hdf : 0
       fast5_raw : 1
       fast5_reads_per_folder : 4000
       fastq_enabled : 0
       fastq_reads_per_file : 4000
       filename : desktop_f51vn4i_20170830_31_mn17984_sequencing_run_colony_a1_rt_36495
       flowcell_type : flo-min107
       local_basecalling : 0
       local_bc_comp_model :
       local_bc_temp_model : template_r9.5_450bps_5mer.jsn
       sample_frequency : 6024
       sequencing_kit : sqk-lsk108
       user_filename_input : colony_a1_rt
     tracking_id
       asic_id : 4176874058
       asic_id_eeprom : 1937759
       asic_temp : 28.392740
       auto_update : 1
       auto_update_source : https://mirror.oxfordnanoportal.com/software/MinKNOW/
       bream_core_version : 1.7.14.1
       bream_is_standard : 0
       bream_ont_version : 1.7.14.1
       device_id : MN17984
       exp_script_name : 71
       exp_script_purpose : sequencing_run
       exp_start_time : 2017-08-30T18:32:48Z
       flow_cell_id : 31
       heatsink_temp : 34.324219
       hostname : DESKTOP-F51VN4I
       installation_type : map
       local_firmware_file : 0
       operating_system : Windows 6.2
       protocol_run_id : 7f6d8043-0051-4b21-b130-5a9ef546ccb7
       protocols_version : 1.7.14
       run_id : 9ae63f7a99bf7f0c1ef6be8ca95e17e979a91d06
       sample_id : colony_a1_rt
       usb_config : 1.1.1_ONT#MinION_fpga_1.1.0#ctrl#Auto
       version : 1.7.14
   Raw
     Channel_307
       Signal
         Data: int16 (34894336,)
       Meta
         description : Grouper
         range : 2131.14
         digitisation : 8192.0
         offset : 13.0
         sample_rate : 6024.0
     ...
   Meta
     duration_samples : 34894336
     sample_rate : 6024
     User
       analysis_conf
         Data: |S3713593 (1,)
"""

class BulkFast5Server:
  def __init__(self, bulk, port):
    self.f = h5py.File(bulk, 'r')
    for r in self.f["IntermediateData"]["Channel_243"]["Reads"]:
        if int(r[1]) == 86:
            print(r)
        elif int(r[1]) > 86:
            break
    self.meta = json.loads(self.f["Meta"]["User"]["analysis_conf"][0].decode('utf-8'))
    print(list(self.meta.keys()))
    print(self.meta["channel_states"])
    print()
    print(self.meta["event_detection"])
    print()
    print(self.meta["experiment"])
    print()
    print(self.meta["histograms"])
    print()
    print(self.meta["read_classification"])
    print()
    print(self.meta["read_detection"])

    print("Websocket server started on port {}.".format(port))
    self.ws_server = websockets.serve(self.ws_loop, 'localhost', port)
    asyncio.get_event_loop().run_until_complete(self.ws_server)
    asyncio.get_event_loop().run_forever()


  @asyncio.coroutine
  def ws_loop(self, websocket, path):
    while True:
      cmd = yield from websocket.recv()
      cmd = cmd.lower().strip().split()
      if cmd[0] == "meta":
        data = self.get_meta()
      elif cmd[0] == "signal":
        channel = int(cmd[1])
        st = int(float(cmd[2]))
        en = int(float(cmd[3]))
        data = self.get_signal(channel, st, en)
      elif cmd[0] == "events":
        channel = int(cmd[1])
        st = int(float(cmd[2]))
        en = int(float(cmd[3]))
        data = self.get_events(channel, st, en)
      elif cmd[0] == "reads":
        channel = int(cmd[1])
        st = int(float(cmd[2]))
        en = int(float(cmd[3]))
        data = self.get_reads(channel, st, en)
      else:
        data = None
      yield from websocket.send(json.dumps(data))


  def get_meta(self):
    try:
      data = {
        "duration_samples": int(self.f["Meta"].attrs["duration_samples"]),
        "sample_rate": int(self.f["Meta"].attrs["sample_rate"]),
        "run_id": self.f["UniqueGlobalKey"]["tracking_id"].attrs["run_id"].decode('utf-8'),
        "sample_id": self.f["UniqueGlobalKey"]["tracking_id"].attrs["sample_id"].decode('utf-8')
      }
    except Exception as e:
      return {"error": "{}\n{}".format(e, traceback.format_exc())}
    return data


  def get_signal(self, channel, st, en):
    if en - st > 100000:
      return {"error": "Too large a range requested, keep it under 100k for now"}
    try:
      channel_range = self.f["IntermediateData"]["Channel_{}".format(channel)]["Meta"].attrs["range"]
      channel_digitisation = self.f["IntermediateData"]["Channel_{}".format(channel)]["Meta"].attrs["digitisation"]
      channel_offset = self.f["IntermediateData"]["Channel_{}".format(channel)]["Meta"].attrs["offset"]
      # the events are scaled as though (measured value = (pA / range * digitisation) - offset)
      # so, we reverse this process below to get pA for the signal, which will then match the event values

      state_data = self.f["StateData"]["Channel_{}".format(channel)]["States"]
      state_data = search_subset(state_data, st, en, 0, [2]) # 0 is the index in the struct of the raw index to match
      mux_data = self.f["MultiplexData"]["Channel_{}".format(channel)]["Multiplex"]
      mux_data = search_subset(mux_data, st, en, 0, [1])
      data = {
        "signal": [int((a + channel_offset) / channel_digitisation * channel_range) for a in self.f["Raw"]["Channel_{}".format(channel)]["Signal"][st:en]],
        "state": state_data,
        "mux": mux_data,
        "start": st
      }
    except Exception as e:
      return {"error": "{}\n{}".format(e, traceback.format_exc())}
    return data


  def get_events(self, channel, st, en):
    if en - st > 100000:
      return {"error": "Too large a range requested, keep it under 100k for now"}
    try:
      event_data = self.f["IntermediateData"]["Channel_{}".format(channel)]["Events"]
      event_data = search_subset(event_data, st, en, 0, [1,3,4]) # 0 is the index in the struct of the raw index to match
      data = {
        "events": event_data
      }
    except Exception as e:
      return {"error": "{}\n{}".format(e, traceback.format_exc())}
    return data


  def get_reads(self, channel, st, en):
    if en - st > 100000:
      return {"error": "Too large a range requested, keep it under 100k for now"}
    try:
      read_data = self.f["IntermediateData"]["Channel_{}".format(channel)]["Reads"]
      read_data = search_subset(read_data, st, en, 2, [3,0,1,4,5]) # length, name, read#, event_start, event_end
      data = {
        "reads": read_data
      }
    except Exception as e:
      return {"error": "{}\n{}".format(e, traceback.format_exc())}
    return data


def search_subset(data, st, en, raw_idx, data_idx):
  hi = len(data) - 1
  lo = 0
  while hi - lo > 1:
    mid = lo + int((hi - lo) / 2)
    if data[mid][raw_idx] < st:
      lo = mid
    elif data[mid][raw_idx] > st:
      hi = mid - 1
    else:
      break
  st_idx = lo if data[hi][raw_idx] > st else hi
  en_idx = st_idx
  while en_idx < len(data)-1 and data[en_idx][raw_idx] < en:
    en_idx += 1
  return [[int(d[raw_idx])] + [int(d[i]) if 'int' in str(type(d[i])) else (float(d[i]) if 'float' in str(type(d[i])) else d[i].decode('utf-8')) for i in data_idx] for d in data[st_idx:en_idx+1]]


def main(bulk, port):
  server = BulkFast5Server(bulk, port)


if __name__ == "__main__":
  parser = argparse.ArgumentParser("Serve access to a bulk FAST5 file over a websocket")
  parser.add_argument("bulk", help="Bulk FAST5 file")
  parser.add_argument("--port", help="Websocket port (default: 8080)", type=int, default=8080)
  args = parser.parse_args()
  main(args.bulk, args.port)
