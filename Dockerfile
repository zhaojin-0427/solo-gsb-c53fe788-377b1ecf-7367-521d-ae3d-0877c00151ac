FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app ./app

EXPOSE 8000

# 单进程：各请求内部 CP-SAT 已并行；避免多 worker 互相争抢 CPU。
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
