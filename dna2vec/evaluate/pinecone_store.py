"""
Pinecone store for DNA2Vec embeddings

This usually crashes in the first run but succeeds in the second 0_o
Issues with API timeout on DB creation
Unstable: https://github.com/pinecone-io/pinecone-python-client/issues
"""
import string
import pinecone
import torch
from tqdm import tqdm
import random
from inference_models import EvalModel, Baseline
from typing import Optional
import concurrent.futures


class PineconeStore:
    def __init__(
        self,
        device,
        index_name: str,
        metric: str = "cosine",
        model_params = None,
        baseline: bool = False,
        baseline_name: Optional[str] = None
    ):
        if model_params is None and not baseline:
            raise ValueError("Model params are empty.")
        if baseline:
            self.model = Baseline(
                option = baseline_name,
                device = device
            )
        else:
            self.model = EvalModel(
                model_params["tokenizer"],
                model_params["model"],
                model_params["pooling"],
                device = device
            )

        if "config-" in index_name: # premium account
            self.api_key = "ded0a046-d0fe-4f8a-b45c-1d6274ad555e"
            self.environment = "us-west4-gcp"
            
        else:
            raise NotImplementedError("Name not identified.")
        
        self.initialize_pinecone_upsertion(metric, index_name)
        self.index_name = index_name

    def initialize_pinecone_upsertion(
        self, 
        metric: str, 
        index_name: str
    ):
        pinecone.init(
            api_key=self.api_key, 
            environment=self.environment
        )

        # only create index if it doesn't exist
        if index_name not in pinecone.list_indexes():
            print(f"Creating new index, {index_name}")
            
            try:
                dimension = self.model.get_sentence_embedding_dimension()
            except:
                dimension = 384
                
            pinecone.create_index(
                name=index_name,
                dimension=dimension,
                metric=metric,
                pod_type="s1.x4"
            )

        # now connect to the index
        self.index = pinecone.GRPCIndex(index_name)

    @staticmethod
    def batched_data_generator(file_path, batch_size):
        """Generator to stream many snippets from text file

        Args:
            file_path (str): Location of file to stream individual snippets
            batch_size (int): stream load

        Yields:
            list(str): list of strings
        """
        import pickle
        
        namespace = ""
        with open(file_path, "rb") as f:
            
            list_of_objects = pickle.load(f)
            batch = []
            
            for unit in list_of_objects:
                
                if namespace == "":
                    namespace = unit["metadata"]
                    
                if namespace != unit["metadata"] or len(batch) >= batch_size:
                    yield batch, namespace
                    batch = [unit]
                    namespace = unit["metadata"]
                else:
                    batch.append(unit)

            if batch:
                yield batch, namespace


    @staticmethod
    def generate_random_string(length: str = 20):
        letters = string.ascii_lowercase
        return "".join(random.choice(letters) for _ in range(length))



    def trigger_pinecone_upsertion(self, file_paths: list, 
                                   batch_size: int = 100, add_namespace=False):
        from tqdm import tqdm
        
        for file_path in file_paths:
            batches = PineconeStore.batched_data_generator(file_path, batch_size)

            for _, (batch,namespace) in tqdm(enumerate(batches)):
                ids = [PineconeStore.generate_random_string() for _ in range(len(batch))]

                # create metadata batch - we can add context here
                metadatas = batch
                texts = [text["text"] for text in batch]
                # create embeddings
                xc = self.model.encode(texts)

                # create records list for upsert
                records = zip(ids, xc, metadatas)
                # upsert to Pinecone
                if add_namespace:
                    self.index.upsert(vectors=records, namespace=namespace)
                else:
                    self.index.upsert(vectors=records)

        # check number of records in the index
        self.index.describe_index_stats()



    # def query(self, query, top_k=5):  # consider batching if slow
    #     # create the query vector
    #     xq = self.model.encode(query).tolist()
    #     # now query
    #     xc = self.index.query(xq, top_k=top_k, include_metadata=True)
    #     return xc

    
    def query_batch(self, queries, indices, top_k=5, hotstart_list=None, meta_dict=None, prioritize=False):  # consider batching if slow
        
        # create the query vector
        xqs = self.model.encode(queries).tolist()
        all_results = []
        
        def query_single(xq, query, index, single_hotstart):
            if not prioritize:
                xc = self.index.query(xq, top_k=top_k, include_metadata=True)
            else:
                xc = self.index.query(xq, top_k=top_k, include_metadata=True,
                                      namespace = single_hotstart)
            
            xc["query"] = query
            xc["index"] = index
            return xc

        if hotstart_list is None:
            with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
                futures = [executor.submit(query_single, xq, query, index, None) \
                    for xq, query, index in zip(xqs, queries, indices)]
            
            for future in concurrent.futures.as_completed(futures):
                all_results.append(future.result())
                
        else:
            # for xq, query, index, single_hotstart in zip(xqs, queries, indices, hotstart_list):
            #     all_results.append(query_single(xq, query, index, single_hotstart))
            with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
                futures = [executor.submit(query_single, xq, query, index, single_hotstart) \
                    for xq, query, index, single_hotstart in zip(xqs, queries, indices, hotstart_list)]

            for future in concurrent.futures.as_completed(futures):
                all_results.append(future.result())

        return all_results
    
    
    def drop_table(self):  # times out for large data!
        pinecone.delete_index(self.index_name)






if __name__ == "__main__":
    import argparse
    import random

    # fmt: off
    parser = argparse.ArgumentParser(description="Pinecone parse and upload.")
    parser.add_argument('--reupload', help="Should we reupload the data?", type=str, choices=['y','n'])
    parser.add_argument('--drop', help="Should we drop the data?", type=str, choices=['y','n'])
    parser.add_argument('--inputpath', help="Splice input file", type=str)
    parser.add_argument('--indexname', help="Index name", type=str)
    
    args = parser.parse_args()
    # fmt: on

    random.seed(42)

    pinecone_obj = PineconeStore(
        device="cuda:3", 
        index_name=args.indexname
    )

    # if args.reupload == "y":
    #     pinecone_obj.trigger_pinecone_upsertion(
    #         file_path=args.inputpath
    #     )

    if args.drop == "y":
        pinecone_obj.drop_table()
